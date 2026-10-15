
## contains the functions for creating, printing, and saving results of simulation

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import cairosvg
from graphviz import Digraph
from sympy import Matrix
import os
import re
import io

class RESULTS():

	def __init__(self, reactor, reactor_type, qss, reaction_type):
		self.reactor = reactor
		self.reactor_type = reactor_type
		self.qss = qss
		self.reaction_type = reaction_type
		self.scheme = reactor.scheme

	def elementary_equations(self):
		if self.reaction_type == "net":
			lines = ["\nPROVIDED NET EQUATIONS:"]
			for reaction in self.scheme._reactions:
				lines.append(f"{reaction['name'] + ':':<6}{reaction['equation']}")

			net_equations_data = "\n".join(lines)
			return net_equations_data

		elif self.reaction_type == "elementary":
			lines = ["\nELEMENTARY EQUATIONS:"]
			for reaction in self.scheme._reactions:
				lines.append(f"{reaction['name'] + ':':<6}{reaction['equation']}")

			elementary_equations_data = "\n".join(lines)
			return elementary_equations_data

	def overall_equations(self):
		# exceptions to printing
		if len(self.scheme._labels_intermediates) == 0:
			return "\nOVERALL REACTIONS: not available as there are no reactive intermediate species"
		if len(self.scheme._reactions) == 1:
			return "\nOVERALL REACTIONS: not available as there is only one reaction"

		reaction_number = 1
		lines = ["\nOVERALL REACTIONS:"]

		# matrix with only the reactive int stoichiometry
		ri_indices = [list(self.scheme._labels).index(label) for label in self.scheme._labels_intermediates]
		SM_ri = Matrix(self.scheme._S).T[:, ri_indices] 

		# finds combinations of reactions that cancel reactive ints
		vectors = SM_ri.T.nullspace() 

		# reproduces combinations as overall equations
		for vec in vectors:
			reaction = (vec.T * Matrix(self.scheme._S).T)[0, :]
			reactants = []
			products = []

			for coeff, species in zip(reaction, list(self.scheme._labels)):
				if coeff < 0:
					if abs(coeff) == 1:
						reactants.append(species)
					else:
						reactants.append(f"{float(abs(coeff)):g} {species}")
				elif coeff > 0:
					if coeff == 1:
						products.append(species)
					else:
						products.append(f"{float(coeff):g} {species}")

			# disregards reactions with no reactants or products
			if not reactants or not products:
				continue 

			

			lines.append(f"{str(reaction_number) + ':':<6}{' + '.join(reactants)} -> {' + '.join(products)}")
			reaction_number += 1

		overall_equation_data = "\n".join(lines)
		return overall_equation_data

	def raw_data(self):
		return self.reactor.result

	def final_concentrations(self):
		units = "mM"
		lines = [f"\nFINAL STABLE CONCENTRATIONS ({units}):"]

		for species in list(self.scheme._labels_stable) + list(self.scheme._labels_buffered):
			final_conc = self.reactor.result[species]["y"][-1]
			if final_conc < 1e-10:
				final_conc = 0
			lines.append(f"{species + ':':<8}{final_conc:g} {units}")

		final_concs = "\n".join(lines)
		return final_concs

	def plot_seperate_results(self, species=None, x_limits=None, y_limits=None, filename=None, show=True):
		for key, i in self.reactor.result.items():
			if key in species:
				plt.figure(figsize=(9, 5))

				plt.plot(i['x'], i['y'], lw=2, label=key)

				plt.title(f"Concentration Profile: {key}", fontsize=18)
				plt.xlabel(r'$t$', fontsize=16)
				plt.ylabel(r'$c_i$', fontsize=16)

				if x_limits:
					plt.xlim(x_limits[0],x_limits[1])

				if y_limits is None:
					max_val = np.max(i['y'])
					min_val = np.min(i['y'])
					y_limits=[min_val,max_val]

				plt.ylim(y_limits[0],y_limits[1])

				plt.legend(frameon=False, fontsize=14)

				plt.tight_layout()

				if filename:
					# Add the species name to the filename
					base, ext = os.path.splitext(filename)
					species_filename = f"{base}_{key}{ext}"
					plt.savefig(species_filename, dpi=300, bbox_inches="tight")

				if show:
					plt.show()

	def plot_results(self, filename=None, show=True):
		plt.figure(figsize=(9, 5))
		plt.suptitle("Concentration Profiles", fontsize=18)

		for key, i in self.reactor.result.items():
			plt.subplot(1,2,1)
			line = plt.plot(i['x'], i['y'], lw=2, label=key)

			plt.subplot(1,2,2)
			if key in self.scheme._labels_intermediates:
				plt.plot(i['x'], i['y'], lw=2, label=key,color=line[0].get_color())

		plt.subplot(1,2,1)
		plt.title("All Species", fontsize=16)
		#plt.yscale("log")
		plt.ylim(-5, 70)
		plt.legend(frameon=False, fontsize=16, bbox_to_anchor=(1.05, 1))
		plt.xlabel(r'$t$', fontsize=16)
		plt.ylabel(r'$c_i$', fontsize=16)

		plt.subplot(1,2,2)
		plt.title("Reactive intermediates", fontsize=16)
		#plt.ylim(-5, 70)
		plt.legend(frameon=False, fontsize=16, bbox_to_anchor=(1.05, 1))
		plt.xlabel(r'$t$', fontsize=16)
		plt.ylabel(r'$c_i$', fontsize=16)

		plt.tight_layout()

		if filename:
			plt.savefig(filename, dpi=300, bbox_inches="tight")
			plt.close()
	
		if show:
			plt.show()

	def plot_scheme(self, filename=None, show=True):
		
		def _parse_species_side(side):
			result = []
			for item in side.split("+"):
				item = item.strip()
				item = re.sub(r"^\d+", "", item)
				result.append(item.strip())
			return result

		def new_repeat_node(species):
			if species not in repeat_counter:
				repeat_counter[species] = 1
			else:
				repeat_counter[species] += 1
			return (species + "_regen_" + str(repeat_counter[species]))

		def follow_species(species):
			if species in visited_species:
				new_id = new_repeat_node(species)
				nodes[new_id] = "[" + species + "]"
				return new_id

			visited_species.add(species)
			nodes[species] = species
			consuming = []

			for reaction in reactions:
				if species in reaction["reactants"]:
					consuming.append(reaction)

			if not consuming:
				return species

			for reaction in consuming:
				buffers = []

				for reactant in reaction["reactants"]:
					if reactant != species and species_type.get(reactant) == "buffered":
						buffers.append(reactant)

				label = ",".join(buffers)

				for product in reaction["products"]:
					if product == species:
						continue
					if species_type.get(product) == "buffered":
						continue
					if species_type.get(product) == "stable" or product not in reacting_species:
						nodes[product] = product
						edges.append((species, product, label))
					else:
						child = follow_species(product)
						edges.append((species, child, label))

			return species

		species_type = {}
		reactions = []
		all_products = set()
		reacting_species = set()
		starts = set()
		
		for s in self.scheme._species:
			species_type[s["name"]] = s["type"]

		for reaction in self.scheme._reactions:
			lhs, rhs = reaction["equation"].split("->")
			reactants = _parse_species_side(lhs.strip())
			products = _parse_species_side(rhs.strip())
			reactions.append({"id": reaction["name"], "reactants": reactants, "products": products})

		for reaction in reactions:
			all_products.update(reaction["products"])
			reacting_species.update(reaction["reactants"])

		for reaction in reactions:
			for species in reaction["reactants"]:
				if species not in all_products and species_type.get(species) != "buffered":
					starts.add(species)
		# starts.add("bu2o2")
		if not starts:
			for reaction in reactions:
				for species in reaction["reactants"]:
					for s in self.scheme._species:
						if s["name"] == species and s["initial concentration"] > 0:
							starts.add(species)

		nodes = {}
		edges = []
		visited_species = set()
		repeat_counter = {}

		for start in starts:
			follow_species(start)

		dot = Digraph("reaction_mechanism", format="svg")
		dot.attr(rankdir="LR", splines="polyline", nodesep="0.6", ranksep="1.2")

		for node, label in nodes.items():
			dot.node(node, label, shape="plaintext")

		for start, end, label in edges:
			if label:
				dot.edge(start, end, label=label)
			else:
				dot.edge(start, end)

		# saves data 
		# if filename = none, shows scheme and saves no data
		svg_bytes = dot.pipe(format="svg")
		# if filename is not none, saves as file
		if filename is not None:
			with open(filename + ".svg", "wb") as f:
				f.write(svg_bytes)

		if show:
			png_bytes = cairosvg.svg2png(bytestring=svg_bytes, dpi=300)
			plt.figure(figsize=(9,5))
			plt.imshow(mpimg.imread(io.BytesIO(png_bytes), format="png"))
			plt.axis("off")
			plt.title("Reaction Mechanism", fontsize=14)
			plt.show()

	def write_info_file(self, info_file_path, overall_runtime, show=False):
		with open(info_file_path, "w") as f:
			f.write("-------------------------------------------------------------\nINPUT DATA:\n-------------------------------------------------------------\n")
			f.write("\nREACTIONS:\n")
			for line in self.scheme._reaction_lines:
				f.write(line + "\n")
			f.write("\nSPECIES:\n")
			for line in self.scheme._species_lines:
				f.write(line + "\n")
			f.write(f"\nREACTOR TYPE: {self.reactor_type}")
			f.write(f"\nQSS: {self.qss}")
			if self.reactor_type == "batch" or self.reactor_type == "cstr":
				f.write(f"\nREACTION TIME: {self.scheme._reaction_time} s")
			if self.reactor_type == "pfr" or self.reactor_type == "cstr":
				f.write(f"\nRESIDENCE TIME: {self.scheme._residence_time} s")
			f.write("\n\n-------------------------------------------------------------\nOUTPUT:\n-------------------------------------------------------------\n")
			f.write(self.scheme._stoichiometric_matrix() + "\n")
			f.write(self.overall_equations() + "\n")
			f.write(f"\nCALCULATION TIME: {overall_runtime:.3f} s\n")
			f.write(self.final_concentrations() + "\n")
			f.write("\n\n\n")

		if show == True:
			os.startfile(info_file_path)

	def export_results(self, file_name, save_plot=True, save_raw_data=True, save_info=True, save_scheme=True):
		if save_plot == False and save_raw_data == False and save_info == False and save_scheme == False:
			return
		else:
			base_directory = f"RS_results/RS_results_{file_name}"
			directory_path = base_directory
			dict_count = 1

			while os.path.exists(directory_path):
				directory_path = f"{base_directory}_{dict_count}"
				dict_count += 1

			os.makedirs(directory_path)

			plot_file_path = os.path.join(directory_path, f"RS_plot_{file_name}.png")
			info_file_path = os.path.join(directory_path, f"RS_info_{file_name}.txt")
			raw_data_file_path = os.path.join(directory_path, f"RS_raw_data_{file_name}.txt")
			scheme_svg_path = os.path.join(directory_path, f"RS_scheme_{file_name}")
			scheme_png_path = os.path.join(directory_path, f"RS_scheme_{file_name}.png")

			if save_plot == True:
				self.plot_results(filename=plot_file_path, show=False)

			if save_info == True:
				self.write_info_file(info_file_path, self.reactor.overall_runtime)

			if save_raw_data == True:
				with open(raw_data_file_path, "w") as f:
					f.write("-------------------------------------------------------------\nRESULTS DATA:\n-------------------------------------------------------------\n\n")
					f.write(str(self.raw_data()))
					
			if save_scheme:
				self.plot_scheme(filename=scheme_svg_path, show=False)
				cairosvg.svg2png(url=scheme_svg_path + ".svg", write_to=scheme_png_path, dpi=600)

	def print_info(self):
		print(f"\nREACTOR TYPE: {self.reactor_type}")
		print(f"QSS: {self.qss}")
		print(f"REACTION TYPE: {self.reaction_type}")
		print(self.elementary_equations())
		print(self.overall_equations())
		print(f"\nCALCULATION TIME: {self.reactor.overall_runtime:.3f} s")
		print(self.final_concentrations())
