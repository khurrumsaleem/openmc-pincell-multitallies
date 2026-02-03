from math import pi
import openmc
import openmc.deplete
import matplotlib.pyplot as plt
import numpy as np
import os
temperature_fuel = 626.85 + 273.15         # 900 Kelvin
temperature_uni = 326.85 + 273.15          # 600 Kelvin

# Cross-section path
os.environ["OPENMC_CROSS_SECTIONS"] = "/home/emil/openmc/Cross_Section_Libraries/endfb-viii.0-hdf5/cross_sections.xml"
###############################################################################
#                              Define materials
###############################################################################

uo2 = openmc.Material(name='UO2 fuel at 2.4% wt enrichment', temperature=temperature_fuel)
uo2.set_density('g/cm3', 10.29769)
uo2.add_element('U', 1., enrichment=2.4)
uo2.add_element('O', 2.)


helium = openmc.Material(name='Helium for gap',temperature=temperature_uni)
helium.set_density('g/cm3', 0.001598)
helium.add_element('He', 2.4044e-4)

zircaloy = openmc.Material(name='Zircaloy 4',temperature=temperature_uni)
zircaloy.set_density('g/cm3', 6.55)
zircaloy.add_element('Sn', 0.014, 'wo')
zircaloy.add_element('Fe', 0.00165, 'wo')
zircaloy.add_element('Cr', 0.001, 'wo')
zircaloy.add_element('Zr', 0.98335, 'wo')

borated_water = openmc.Material(name='Borated water',temperature=temperature_uni)
borated_water.set_density('g/cm3', 0.740582)
borated_water.add_element('B', 4.0e-5)
borated_water.add_element('H', 5.0e-2)
borated_water.add_element('O', 2.4e-2)
borated_water.add_s_alpha_beta('c_H_in_H2O')

###############################################################################
# 3. GEOMETRY
###############################################################################
z_bot = openmc.ZPlane(z0=-150)
z_top = openmc.ZPlane(z0=150)
fuel_cyl = openmc.model.RightCircularCylinder(
    center_base=(0.0, 0.0, -150),
    height=300,
    radius=0.47
)

gap_cyl = openmc.model.RightCircularCylinder(
    center_base=(0.0, 0.0, -150),
    height=300,
    radius=0.48
)

clad_cyl = openmc.model.RightCircularCylinder(
    center_base=(0.0, 0.0, -150),
    height=300,
    radius=0.54
)
box = openmc.model.RectangularParallelepiped(
    xmin=-0.75, xmax=0.75,
    ymin=-0.75, ymax=0.75,
    zmin=-160, zmax=160,
    boundary_type='reflective'  # Serpent 
)
# --- Pin cells
fuel_cell = openmc.Cell(
    region= -fuel_cyl & +z_bot & -z_top,
    fill=uo2
)
gap_cell = openmc.Cell(
    region= +fuel_cyl & -gap_cyl & +z_bot & -z_top,
    fill=helium
)
clad_cell = openmc.Cell(
    region= +gap_cyl & -clad_cyl & +z_bot & -z_top,
    fill=zircaloy
)
coolant_cell = openmc.Cell(
    region= +clad_cyl & -box,
    fill=borated_water
)
bottom_reflector = openmc.Cell(
    region= -box & -z_bot,
    fill=borated_water
)

top_reflector = openmc.Cell(
    region= -box & +z_top,
    fill=borated_water
)
outside = openmc.Cell(
    region= +box,
    fill=borated_water
)
geometry = openmc.Geometry([fuel_cell,gap_cell,clad_cell,coolant_cell,bottom_reflector,top_reflector,outside])
height=300
fuel_radius = 0.47  # cm
uo2.volume = pi * height *fuel_radius**2
###############################################################################
#                     Transport calculation settings
###############################################################################

settings = openmc.Settings()
#settings.batches = 100
#settings.inactive = 50
#settings.particles = 1000
settings.particles = 10000
settings.batches = 500
settings.inactive = 200
settings.temperature = {
    'default': temperature_uni,
    'method': 'interpolation',
}
settings.export_to_xml()
#bounds = [-0.45, -0.45, -140, 0.45, 0.45, 140]
#uniform_dist = openmc.stats.Box(*bounds)
#settings.source = openmc.IndependentSource(
#    space=uniform_dist, constraints={'fissionable': True})
# Create an initial uniform spatial source distribution over fissionable zones
settings.source = openmc.source.Source(space=openmc.stats.Point())
entropy_mesh = openmc.RegularMesh()
entropy_mesh.lower_left = [-0.6, -0.6, -150.0]
entropy_mesh.upper_right = [ 0.6,  0.6,  150.0]
entropy_mesh.dimension = [10, 10, 10]
settings.entropy_mesh = entropy_mesh

###############################################################################
#                              Mesh tallies
###############################################################################

mesh = openmc.RegularMesh()
mesh.dimension   = [100, 100, 100]
mesh.lower_left  = [-0.63, -0.63, -150]
mesh.upper_right = [ 0.63,  0.63,  150]
mesh_filter = openmc.MeshFilter(mesh)

# Define tallies
flux_tally = openmc.Tally(name='flux mesh tally')
flux_tally.filters = [mesh_filter]
flux_tally.scores = ['flux']

fission_tally = openmc.Tally(name='fission mesh tally')
fission_tally.filters = [mesh_filter]
fission_tally.scores = ['fission']

tallies = openmc.Tallies([flux_tally, fission_tally])

###############################################################################
#                   Initialize and run depletion calculation
###############################################################################

model = openmc.Model(geometry=geometry, settings=settings, tallies=tallies)
chain_file = 'chain_simple.xml'
op = openmc.deplete.CoupledOperator(model, chain_file)

#time_steps = [1.0, 1.0, 1.0, 1.0, 1.0]  # days
# cumulative_days = np.linspace(0.0, 360.0, 19)
# time_steps = np.diff(cumulative_days)
#burnup_cum = np.array([
#    0.01, 0.25, 0.5, 1.0, 2.0,
#    3.0, 5.0, 7.5, 10.0,
#    12.5, 15.0, 17.5, 20.0
#])
burnup_cum = np.linspace(0.01, 20.0, 50)
#burnup_cum = np.linspace(0.1, 50.0, 46)
burnup = np.diff(burnup_cum, prepend=0.0)
power = 60000.0  # W/cm (for 2D simulation)
integrator = openmc.deplete.PredictorIntegrator(op, burnup, power, timestep_units='MWd/kg')
print("Initial heavy metal mass [g] =", op.heavy_metal)
print("Initial heavy metal mass [kg] =", op.heavy_metal / 1000)

#integrator.integrate()

