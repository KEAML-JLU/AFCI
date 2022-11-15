#!/usr/bin/env python

# This module helps creating Gantt from a dictionary or a text file.
# Output formats are a Matplotlib chart or a LaTeX code (using pgfgantt).

import random
import numpy as np
import pdb
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from matplotlib import colors as mcolors
import matplotlib.patches as mpatches
colors = []

for name, hex in mcolors.cnames.items():
    colors.append(name)


def parse_data(file):
    try:
        textlist = open(file).readlines()
    except:
        return

    data = {}

    for tx in textlist:
        if not tx.startswith('#'):
            splitted_line = tx.split(',')
            machine = splitted_line[0]
            operations = []

            for op in splitted_line[1::]:
                label = op.split(':')[0].strip()
                l = op.split(':')[1].strip().split('-')
                start = int(l[0])
                end = int(l[1])
                operations.append([start, end, label])

            data[machine] = operations
    return data

def colorcode(jobnum: int) -> str:
    """
    负责给甘特图上色
    """
    if jobnum == 0:
        return 'rosybrown'
    elif jobnum == 1:
        return 'indianred'
    elif jobnum == 2:
        return 'tomato'
    elif jobnum == 3:
        return 'lightsalmon'
    elif jobnum == 4:
        return 'peru'
    elif jobnum == 5:
        return 'darkorange'
    elif jobnum == 6:
        return 'tan'
    elif jobnum == 7:
        return 'gold'
    elif jobnum == 8:
        return 'olive'
    elif jobnum == 9:
        return 'y'
    elif jobnum == 10:
        return 'greenyellow'
    elif jobnum == 11:
        return 'darkseagreen'
    elif jobnum == 12:
        return 'limegreen'
    elif jobnum == 13:
        return 'aquamarine'
    elif jobnum == 14:
        return 'lightseagreen'
    elif jobnum == 15:
        return 'darkcyan'
    elif jobnum == 16:
        return 'deepskyblue'
    elif jobnum == 17:
        return 'deeppink'
    elif jobnum == 18:
        return 'palevioletred'
    elif jobnum == 19:
        return 'pink'
    else:
        return 'black'
def draw_chart(data, svg_name):
    nb_row = len(data.keys())

    pos = np.arange(0.5, nb_row * 0.5 + 0.5, 0.5)

    fig = plt.figure(figsize=(20, 8))
    ax = fig.add_subplot(111)

    index = 0
    max_len = []

    for machine, operations in sorted(data.items()):
        for op in operations:
            max_len.append(op[1])
            # c = random.choice(colors)
            temp = op[2].split('-')

            c = colorcode(int(temp[0]))
            rect = ax.barh((index * 0.5) + 0.5, op[1] - op[0], left=op[0], height=0.3, align='center',
                           edgecolor=c, color=c, alpha=0.8)

            # adding label
            width = int(rect[0].get_width())
            # Str = "J{}".format(op[2])
            Str = temp[1]
            xloc = op[0] + 0.50 * width
            clr = 'black'
            align = 'center'

            yloc = rect[0].get_y() + rect[0].get_height() / 2.0
            ax.text(xloc, yloc, Str, horizontalalignment=align,
                            verticalalignment='center', color=clr, weight='bold',
                            clip_on=True)
        index += 1
    # pdb.set_trace()
    ax.set_ylim(ymin=-0.1, ymax=nb_row * 0.5 + 0.5)
    ax.grid(color='gray', linestyle=':')
    ax.set_xlim(0, max(10, max(max_len)))

    labelsx = ax.get_xticklabels()
    plt.setp(labelsx, rotation=0, fontsize=10)

    locsy, labelsy = plt.yticks(pos, data.keys())
    plt.setp(labelsy, fontsize=14)

    font = font_manager.FontProperties(size='small')
    ax.legend(loc=1, prop=font)

    ax.invert_yaxis()
    color = ['rosybrown','indianred','tomato','lightsalmon','peru','darkorange','tan','gold','olive','y','greenyellow','darkseagreen','limegreen','aquamarine','lightseagreen','darkcyan','deepskyblue','deeppink','palevioletred','pink']
    labels = [i for i in range(1,21)]
    patches = [ mpatches.Patch(color=color[i], label='{:s}'.format(str(labels[i])) ) for i in range(len(color)) ]
    ax=plt.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.height , box.width* 0.8])
    ax.legend(handles=patches, loc = 'center left', bbox_to_anchor=(1,0.5), ncol=1)
    plt.title("Job Shop Solution")
    # plt.title("Flexible Job Shop Solution")
    plt.savefig(svg_name + '.svg')
    plt.show()


def export_latex(data):
    max_len = []
    head = """
\\noindent\\resizebox{{\\textwidth}}{{!}}{{
\\begin{{tikzpicture}}[x=.5cm, y=1cm]
\\begin{{ganttchart}}{{1}}{{{}}}
[vgrid, hgrid]{{{}}}
\\gantttitle{{Flexible Job Shop Solution}}{{{}}} \\\\
\\gantttitlelist{{1,...,{}}}{{1}} \\\\
"""
    footer = """
\\end{ganttchart}
\\end{tikzpicture}}\n
    """
    body = ""
    for machine, operations in sorted(data.items()):
        counter = 0
        for op in operations:
            max_len.append(op[1])
            label = "O$_{{{}}}$".format(op[2].replace('-', ''))
            body += "\\Dganttbar{{{}}}{{{}}}{{{}}}{{{}}}".format(machine, label, op[0]+1, op[1])
            if counter == (len(operations) - 1):
                body += "\\\\ \n"
            else:
                body += "\n"
            counter += 1

    lenM = max(10, max(max_len))
    print(head.format(lenM, lenM, lenM, lenM))
    print(body)
    print(footer)


if __name__ == '__main__':
    fname = r"test.txt"
    draw_chart(parse_data(fname))
