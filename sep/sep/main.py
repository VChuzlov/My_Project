import json
from flow import flow
from separator2f import separator2f
from separator3f import separator3f
from compressor import compressor
from energy_flow import energy_flow
from chiller import chiller
from valve import valve

with open(r'.\input_sep.json', 'r') as input_file:
    input = json.load(input_file)

flow1 = flow('flow1')
flow2 = flow('flow2')
flow2 = flow('flow3')
flow3 = flow('flow4')
flow4 = flow('flow5')
flow5 = flow('flow6')
flow6 = flow('flow7')
flow7 = flow('flow8')
flow8 = flow('flow9')
flow9 = flow('flow10')
energy_flow1 = energy_flow('energyflow1')
energy_flow2 = energy_flow('energyflow2')

flow1.pressure = input['flow1_pressure']
flow1.temperature = input['flow1_temperature']
flow1.mass_flow = input['flow1_mass_flow']
flow1.composition = input['flow1_composition']
flow1.calculate()

separator2f1 = separator2f('separator2f1', flow1)
flow2, flow3 = separator2f1.calculate()
flow2.name, flow3.name = 'flow2', 'flow3' #

compressor1 = compressor('compressor1', flow2, 8200)
flow4, energy_flow1.energy = compressor1.calculate()
flow4.name = 'flow4'

chiller1 = chiller("chiller1", flow4, 10)
flow5, energy_flow2.energy = chiller1.calculate()
flow5.name= 'flow5'

valve1 = valve("valve1", flow5, 6200)
flow6 = valve1.calculate()
flow6.name = 'flow6'

separator3f1 = separator3f("separator3f1", flow6)
flow7, flow8, flow9 = separator3f1.calculate()
flow7.name, flow8.name, flow9.name = 'flow7', 'flow8', 'flow9'

results = {}
for object in [flow1, flow2, flow3, flow4, flow5,
                flow6, flow7, flow8, flow9, separator2f1,
                compressor1, chiller1, valve1, separator3f1]:
    results.update(object.results())

with open(r'.\result_sep.json', "w") as result_file:
    result_file.write(json.dumps(results))
pass