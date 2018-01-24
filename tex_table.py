# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 19 14:12:08 2017

Write the result into tex table.

@author: Neo
"""

# ----------- Horizontal Table -------------------


def htable(names, estimates, errors, fout):
    '''save estmates and corresponding formal errors into a horizontal table.
    '''
    nameline = '$%10s$' % names[0]
    # dataline = '%+8.3f $\\pm$ %8.3f' % (estimates[0], errors[0])
    dataline = '$%+8.1f \\pm %8.1f$' % (estimates[0], errors[0])
    for i in range(1, len(names)):
        nameline += '  &$%10s$' % names[i]
        # dataline += '  &%+8.3f $\\pm$ %8.3f' % (estimates[i], errors[i])
        dataline += '  &$%+8.1f \\pm %8.1f$' % (estimates[i], errors[i])
    nameline += ' \\\\'
    dataline += ' \\\\'
    print(nameline + '\n \\hline \n' + dataline, file=fout)
# ----------- END -------------------
# ----------- Vertical Table -------------------


def vtable(names, estimates, errors, fout):
    '''save estmates and corresponding formal errors into a vertical table.
    '''
    for i in range(len(names)):
        # print('$%10s$  &%+8.3f $\\pm$ %8.3f \\\\'\
        print('$%10s$  &$%+8.1f \\pm %8.1f$ \\\\'
              % (names[i], estimates[i], errors[i]), file=fout)
# ----------- END -------------------
# ----------- Correlation coefficients -------------------


def cor_table(names, cor, fout):
    print('## correlation coefficients.', file=fout)
    num = len(names)
# Heading
    headline = ' ' * 10
    for name in names:
        headline += '  &$%10s$' % name
    headline += '\\\\ \\hline'
    print(headline, file=fout)
# Coefficient
    for i, name in enumerate(names):
        headpart = '  $%10s$' % name
        blankpart = '  &      ' * i
        datapart = ''
        for j in range(i, num):
            datapart += '  &$%+5.2f$' % cor[i, j]
        print(headpart + blankpart + datapart + '\\\\', file=fout)


def write_result_deg1(x1name, x1, sig1, corr1, flog):
    htable(x1name, x1, sig1, flog)
    cor_table(x1name, corr1, flog)
# -----------------------------------


def write_result_deg2(x1name, x2name, x2, sig2, corr2, flog):
    x21, sig21 = x2[: 6], sig2[: 6]
    x22, sig22 = x2[6:], sig2[6:]
    htable(x1name, x21, sig21, flog)
    htable(x2name, x22, sig22, flog)
    cor_table(x1name + x2name, corr2, flog)
# ----------- END -------------------
