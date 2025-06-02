# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <olivier.dazel@univ-lemans.fr> and <alan.geslain@univ-lemans.fr> wrote this file.  
# As long as you retain this notice you can do whatever you want with this stuff. 
# If we meet some day, and you think this stuff is worth it, you can buy us a beer in return.       
# Olivier DAZEL and Alan GESLAIN
# ----------------------------------------------------------------------------

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons, CheckButtons
import matplotlib.image as mpimg
import numpy as np
import os

from material import possible_fluid_models, possible_solid_models
from Impedance_tube import Impedance_tube, possible_output_types

frequency = np.linspace(1.156250E+2, 4.323437E+3, 500)

class parameter():
    def __init__(self, symbol, value, isactive=False, min_value=None, max_value=None):
        self.symbol = symbol
        self.value = value
        self.isactive = isactive
        self.min_value = min_value if min_value is not None else 0.1 * value
        self.max_value = max_value if max_value is not None else 2 * value

params = {
    'Thickness': parameter('d [cm]', 5.0, True, 1, 10),
    'Porosity': parameter(r"$\phi$ $[\emptyset]$", 0.99, True, 0.8, 1.0),
    'Resistivity': parameter(r'$\sigma \, [Ns/m4]$',9116, True,1000,50000),
    'Tortuosity': parameter(r'$\alpha_\infty [\emptyset]$',1.00, True, 1.0, 2.0),   
    'Viscous characteristic length': parameter(r'$\Lambda \,[\mu m]$',145, True),
    'Thermal characteristic length': parameter(r'$\Lambda'+"'"+r' [\mu m]$',366, True),
    'Static thermal permeability': parameter(r'k'+"'"+'$_0$'+" [$m^2$]",10e-8, True),
    'Static tortuosity': parameter(r'$\alpha_0$ $[\emptyset]$',4/3, True,1,2),
    'Frame density': parameter(r'$\rho_1$ $[kg/m^3]$',10, True),
    'Young\'s modulus': parameter(r'E [Pa]',57313, True),
    'Poisson\'s ratio': parameter(r'$\nu$ $[\emptyset]$',0.23, True),
    'Loss factor': parameter(r'$\eta$ $[\emptyset]$',0.0978, True),
    'frequency': frequency,
    'termination_condition': 'Rigid',
    'fluid_model': 'Delany Bazley',
    'solid_model': 'eqf'
}

pb = Impedance_tube(params)

fig = plt.figure(figsize=(11, 8))
ax = fig.add_axes([0.3, 0.25, 0.65, 0.4])

# Indicator
ax_indicator = fig.add_axes([0.15, 0.86, 0.18, 0.08])
indicator_radio = RadioButtons(ax_indicator, possible_output_types)
ax_indicator.text(0.5, 1.1, 'Indicator', ha='center', transform=ax_indicator.transAxes,
                fontsize=11, fontweight='bold')

# --- CheckButton "Abs mea" ---
#ax_abs_mea = fig.add_axes([0.03, 0.7, 0.17, 0.05])
ax_abs_mea = fig.add_axes([0.15, 0.825, 0.18, 0.035])
check_abs_mea = CheckButtons(ax_abs_mea, ['Show measurement'], [False])

# Fluid models
ax_model_fluid = fig.add_axes([0.45, 0.76, 0.25, 0.18], facecolor='lightblue')
model_check_fluid = CheckButtons(ax_model_fluid, possible_fluid_models, [False]*len(possible_fluid_models))
ax_model_fluid.text(0.5, 1.05, 'Fluid phase model', ha='center',
                    transform=ax_model_fluid.transAxes, fontsize=11, fontweight='bold', color='blue')

# Solid models
ax_model_solid = fig.add_axes([0.78, 0.84, 0.17, 0.1], facecolor='orange')
model_check_solid = CheckButtons(ax_model_solid, possible_solid_models, [False]*len(possible_solid_models))
ax_model_solid.text(0.5, 1.05, 'Solid frame', ha='center',
                    transform=ax_model_solid.transAxes, fontsize=11, fontweight='bold', color='orange')

# Sliders
sliders = {}
slider_start_y = 0.6
slider_height = 0.02
slider_step = 0.025
slider_width = 0.07
slider_x = 0.1
slider_y = slider_start_y

for name, p in params.items():
    if isinstance(p, parameter) and p.isactive:
        ax_slider = fig.add_axes([slider_x, slider_y, slider_width, slider_height])
        slider = Slider(ax_slider, p.symbol, p.min_value, p.max_value, valinit=p.value)
        slider.label.set_size(9)
        sliders[name] = slider
        slider_y -= slider_step

# Logos
logo_ax = fig.add_axes([0.2, 0.05, 0.1, 0.1])
logo_ax.axis('off')
logo_img = mpimg.imread('logo_irp.png')
logo_ax.imshow(logo_img)

laum_ax = fig.add_axes([0.5, 0.0, 0.3, 0.2])
laum_ax.axis('off')
laum_img = mpimg.imread('logo_laum.png')
laum_ax.imshow(laum_img)


# Initial state
model_check_fluid.set_active(2)
model_check_solid.set_active(0)
indicator_radio.set_radio_props({"color": ["black", "black", "black",]})


# --- Fonction de chargement des données mesurées ---
def load_measurement():
    try:
        # path = os.path.join('measure', 'Absorption', 'Abs_meas.txt')
        with open("measurement.txt", 'r', encoding='utf-8') as f:
            lines = f.readlines()
        # Remplacer les virgules par des points, ignorer les lignes vides ou commentées
        lines = [line.replace(',', '.').strip() for line in lines if line.strip() and not line.startswith('%')]
        data = np.array([list(map(float, line.split())) for line in lines])
        measurement = {'frequency': data[:, 0], 'Absorption': data[:, 1], "R": data[:, 2]+1j*data[:, 3], "Surface Impedance": data[:, 4]+1j*data[:,5]}
        return measurement
    except Exception as e:
        print("Erreur de lecture meas.txt :", e)
        return None

# --- Mise à jour des tracés ---
def update(event):

    selected_fluid_models = [label for label, state in zip(possible_fluid_models, model_check_fluid.get_status()) if state]
    selected_solid_models = [label for label, state in zip(possible_solid_models, model_check_solid.get_status()) if state]
    
    # Colors of the sliders 
    for key in sliders.keys():
        sliders[key].poly.set_fc("red")
    # Some are always green     
    sliders["Thickness"].poly.set_fc("green")
    sliders["Porosity"].poly.set_fc("green")
    sliders["Resistivity"].poly.set_fc("green")

    flag = any(m in ["Johnson-Champoux-Allard", "Johnson-Lafarge", "Pride-Lafarge"] for m in selected_fluid_models)
    if flag:
        sliders["Tortuosity"].poly.set_fc("green")
        sliders["Tortuosity"].poly.set_fc("green")
        sliders["Viscous characteristic length"].poly.set_fc("green")
        sliders["Thermal characteristic length"].poly.set_fc("green")
    flag = any(m in ["Johnson-Lafarge", "Pride-Lafarge"] for m in selected_fluid_models)
    if flag:
        sliders["Static thermal permeability"].poly.set_fc("green")
    if "Pride-Lafarge" in selected_fluid_models:
        sliders["Static tortuosity"].poly.set_fc("green")

    flag = any(m in ["Biot", "Limp"] for m in selected_solid_models)
    if flag:
        sliders["Frame density"].poly.set_fc("green")
    if "Biot" in selected_solid_models:
        sliders["Young\'s modulus"].poly.set_fc("green")
        sliders["Poisson\'s ratio"].poly.set_fc("green")
        sliders["Loss factor"].poly.set_fc("green")

    for name, slider in sliders.items():
        params[name].value = slider.val
    
    params["termination_condition"] = "rigid"

    indicator = indicator_radio.value_selected
    ax.clear()
    

    if indicator == 'Absorption':
        ax.set_title("Absorption Coefficient", fontsize=12, fontweight='bold')
        for fluid_model in selected_fluid_models:
            params['fluid_model'] = fluid_model  # Set the current fluid model in the pb instance
            for solid_model in selected_solid_models:
                params['solid_model'] = solid_model # Set the current solid model in the pb instance
                result = pb.run(frequency, params)
                ax.plot(result['frequency'], result['Absorption'], label=f'{fluid_model} - {solid_model}')
        ax.set_ylabel('Absorption Coefficient')
        ax.set_ylim(0, 1.1)

        # Tracer les données mesurées si activé
        if check_abs_mea.get_status()[0]:
            measurement = load_measurement()
            if measurement['frequency'] is not None:
                ax.scatter(measurement['frequency'], measurement['Absorption'], color='black', s=12, marker='o', label='Measurement')

    elif indicator == 'Surface Impedance':
        ax.set_title("Surface Impedance", fontsize=12, fontweight='bold')
        for fluid_model in selected_fluid_models:
            params['fluid_model'] = fluid_model  # Set the current fluid model in the pb instance
            for solid_model in selected_solid_models:
                params['solid_model'] = solid_model # Set the current solid model in the pb instance
                result = pb.run(frequency, params)
                ax.plot(result['frequency'], np.real(result[indicator]), label=f'{fluid_model} - {solid_model} (Re)')
                ax.plot(result['frequency'], np.imag(result[indicator]), label=f'{fluid_model} - {solid_model} (Im)')
        ax.set_ylabel('Surface Impedance')
        # Tracer les données mesurées si activé
        if check_abs_mea.get_status()[0]:
            measurement = load_measurement()
            if measurement['frequency'] is not None:
                ax.scatter(measurement['frequency'], np.real(measurement[indicator]), color='black', s=12, marker='o', label='Measureament (Re)')
                ax.scatter(measurement['frequency'], np.imag(measurement[indicator]), color='gray', s=12, marker='o', label='Measureament (Im)')



    elif indicator == 'R':
        ax.set_title("Reflection Coefficient", fontsize=12, fontweight='bold')
        for fluid_model in selected_fluid_models:
            params['fluid_model'] = fluid_model  # Set the current fluid model in the pb instance
            for solid_model in selected_solid_models:
                params['solid_model'] = solid_model # Set the current solid model in the pb instance
                result = pb.run(frequency, params)
                ax.plot(result['frequency'], np.real(result[indicator]), label=f'{fluid_model} - {solid_model} (Re)')
                ax.plot(result['frequency'], np.imag(result[indicator]), label=f'{fluid_model} - {solid_model} (Im)')
        if check_abs_mea.get_status()[0]:
            measurement = load_measurement()
            if measurement['frequency'] is not None:
                ax.scatter(measurement['frequency'], np.real(measurement[indicator]), color='black', s=12, marker='o', label='Measureament (Re)')
                ax.scatter(measurement['frequency'], np.imag(measurement[indicator]), color='gray', s=12, marker='o', label='Measureament (Im)')
        ax.set_ylabel('Reflection Coefficient')


    ax.set_xlabel('Frequency [Hz]')
    ax.legend(fontsize=9)
    fig.canvas.draw_idle()

# --- Callbacks ---
def register_callbacks():
    model_check_fluid.on_clicked(update)
    model_check_solid.on_clicked(update)
    indicator_radio.on_clicked(update)
    check_abs_mea.on_clicked(update)
    for cb in sliders.values():
        cb.on_changed(lambda val: update(None))


register_callbacks()
update(None)
plt.show()
