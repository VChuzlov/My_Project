NAMES = (
    'Nitrogen', 'CO2', 'Methane', 'Ethane', 'Propane',
    'i-Butane', 'n-Butane', 'i-Pentane', 'n-Pentane'
    'n-Hexane', 'H2O', 'Methanol', 'C5+*', 'C7+*'
    'C11+*', 'C12+*', 'C19+*', 'C21+*', 'C31+*',
    'C36+*', 'F1*', 'F2*', 'F3*', 'F4*',
)

COMP_COUNT = len(NAMES)

MR = (

)

DENSITY = (

)

COEF_MASS_CP = (

)

TC = [
    220 for _ in range(COMP_COUNT)
]

PC = [
    2.15 for _ in range(COMP_COUNT)
]

OMEGA = [
    0.5 for _ in range(COMP_COUNT)
]

V = (
    0.8 for _ in range(COMP_COUNT)
)

R = 8.314