import collections

import sympy
from brian2.units import *

def parameter_lines(filename):
    with open(filename, 'r') as f:
        extracted = []
        include = False
        for line in f.readlines():
            if line.startswith('# Model parameters'):
                include = True
            elif include and line.startswith('#'*80) and len(extracted) > 1:
                break
            if include and not line.startswith('#') and len(line.strip()):
                extracted.append(line)
        return extracted

def parse_line(line):
    try:
        if '#' in line:
            assignment, comment = line.split('#')
        else:
            assignment = line
            comment = ''
        varname, value = assignment.split('=')
    except ValueError:
        print 'Cannot parse line: {}'.format(line)
        return None, None, None
    return varname.strip(), value.strip(), comment.strip()

def parse_lines(lines):
    parameters = {}
    for line in lines:
        varname, value, comment = parse_line(line)
        parameters[varname] = (value, comment)
    return parameters

def parse_file(filename):
    return parse_lines(parameter_lines(filename))

def merge_parameters(all_parameters):
    all_varnames = {key for params in all_parameters.itervalues() for key in params.keys()}
    common = {}
    different = {}
    for varname in all_varnames:
        filenames = []
        values = []
        comments = []
        for filename, parameters in all_parameters.iteritems():
            if varname in parameters:
                filenames.append(filename)
                value, comment = parameters[varname]
                values.append(value)
                comments.append(comment)
        unique_comments = {c for c in comments if c is not None and len(c)}
        if len({c for c in comments if c is not None and len(c)}) > 1:
            print 'WARNING: Different comments for variable {} ({})'.format(varname, unique_comments)
        if len(set(values)) != 1 or len(values) == 1:
            this_parameter = {}
            for value in set(values):
                filenames_for_value = tuple([f for index, f in enumerate(filenames) if values[index] == value])
                this_parameter[filenames_for_value] = value
            # Special case: all simulations share one value, and only one simulation uses another one
            if len(this_parameter) == 2 and any(len(k) == 1 for k in this_parameter):
                for key in this_parameter.keys():
                    if len(key) > 1:
                        common[varname] = (this_parameter[key], comments[0])
                        del this_parameter[key]
            different[varname] = (this_parameter, comments[0])
        else:
            common[varname] = (values[0], comments[0])
    return common, different

def latex_table_common_parameters(parameters):
    start = r'''
    \section*{Common parameters}
    \noindent
    \begin{tabularx}{\textwidth}{llrX}
    \emph{Symbol} & \emph{Name in code} & \emph{Value} & \emph{Description}\\
    '''
    end = r'''
    \end{tabularx}
    '''
    parameter_lines = []
    for parameter, (value, comment) in sorted(parameters.iteritems()):
        code_variable, latex_variable, value = to_latex(parameter, value)
        parameter_lines.append(r'$%s$ & %s & $%s$ & \detokenize{%s}\\' % (latex_variable, code_variable, value, comment))
    return start + '\n'.join(parameter_lines) + end


def to_latex(parameter, value):
    latex_variable = sympy.latex(sympy.sympify(parameter))
    code_variable = r'\lstinline|{}|'.format(parameter)
    try:
        value = sympy.latex(eval(value))
    except NameError:
        # referring to another name, not a numerical value
        vaule = sympy.latex(sympy.sympify(value))
    return code_variable, latex_variable, value


def latex_individual_parameters(parameters):
    parameters_per_file = collections.defaultdict(list)
    for param, (values, comment) in parameters.iteritems():
        for subset, value in values.iteritems():
            for filename in subset:
                parameters_per_file[filename].append((param, value, comment))
    lines = [r'\section*{Individual parameters}']
    for filename, file_parameters in sorted(parameters_per_file.iteritems()):
        lines.append(r'\subsection*{\texttt{%s}}' % filename.replace('_', r'\_'))
        variables = []
        for param, value, comment in file_parameters:
            code_variable, latex_variable, value = to_latex(param, value)
            variables.append('$%s = %s$ (\detokenize{%s})' % (latex_variable, value, comment))
        lines.append('Parameters: ' + ', '.join(sorted(variables)))
        lines.append('')
    return '\n'.join(lines)


if __name__ == '__main__':
    filenames = ['example_1_COBA.py', 'example_2_gchi_astrocyte.py',
                 'example_3_io_synapse.py', 'example_4_synrel.py',
                 'example_4_rsmean.py', 'example_5_astro_ring.py',
                 'example_6_COBA_with_astro.py']
    all_parameters = {}
    for filename in filenames:
        all_parameters[filename] = parse_file(filename)
    common, different = merge_parameters(all_parameters)
    start = r'''
    \documentclass[DIV=15,fontsize=9pt]{scrartcl}
    \usepackage{tabularx,listings}
    \begin{document}
    '''
    end = r'''
    \end{document}'''
    with open('../text/parameters.tex', 'w') as f:
        f.write(start)
        f.write(latex_table_common_parameters(common))
        f.write(latex_individual_parameters(different))
        f.write(end)
