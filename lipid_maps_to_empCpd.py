import sys
import json

list_emp_Cpds = []
prototype_empCpd = {
    'primary_id': None,
    'primary_db': 'LMSD',
    'name': None,
    'neutral_formula': None,
    'neutral_formula_mass': None,
    'SMILES': None,
    'inchikey': None,
    'other_ids': None,
    'category': None,
    'class': None
}

emp_Cpd = dict(prototype_empCpd)
mode = None
for line in open(sys.argv[1]):
    if mode:
        emp_Cpd[mode] = line.rstrip()
        mode = None
    if line.startswith('> <SMILES>'):
        mode = 'SMILES'
    elif line.startswith('> <LM_ID>'):
        mode = "primary_id"
    elif line.startswith('> <NAME>'):
        mode = 'name'
    elif line.startswith('> <CATEGORY>'):
        mode = 'category'
    elif line.startswith('> <MAIN_CLASS>'):
        mode = 'class'
    elif line.startswith('> <EXACT_MASS>'):
        mode = 'neutral_formula_mass'
    elif line.startswith('> <FORMULA>'):
        mode = 'neutral_formula'
    elif line.startswith('$$$$'):
        mode = None
        list_emp_Cpds.append(emp_Cpd)
        print(emp_Cpd)
        emp_Cpd = dict(prototype_empCpd)
    else:
        mode = None

json.dump(list_emp_Cpds, open(sys.argv[2], 'w+'), indent=4)