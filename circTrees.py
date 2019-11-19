#!/usr/bin/env python3

from Bio import Phylo
import numpy as np
import matplotlib as plt
import pylab

def drawMark(tree, label_func=str, do_show=True, show_confidence=True, 
             axes=None, branch_labels=None, label_colors=None, mark=[],
             markColor='r', markWeight='bold', lineColor='black',
             *args, **kwargs):
    try: 
        import matplotlib.pyplot as plt 
    except ImportError: 
        try: 
            import pylab as plt 
        except ImportError: 
            raise MissingPythonDependencyError( 
                "Install matplotlib or pylab if you want to use draw.") 
 
    import matplotlib.collections as mpcollections 
 
    # Arrays that store lines for the plot of clades 
    horizontal_linecollections = [] 
    vertical_linecollections = [] 
 
    # Options for displaying branch labels / confidence 
    def conf2str(conf): 
        if int(conf) == conf: 
            return str(int(conf)) 
        return str(conf) 
    if not branch_labels: 
        if show_confidence: 
            def format_branch_label(clade): 
                if hasattr(clade, 'confidences'): 
                    # phyloXML supports multiple confidences 
                    return '/'.join(conf2str(cnf.value) 
                                    for cnf in clade.confidences) 
                if clade.confidence is not None: 
                    return conf2str(clade.confidence) 
                return None 
        else: 
            def format_branch_label(clade): 
                return None 
    elif isinstance(branch_labels, dict): 
        def format_branch_label(clade): 
            return branch_labels.get(clade) 
    else: 
        if not callable(branch_labels): 
            raise TypeError("branch_labels must be either a " 
                            "dict or a callable (function)") 
        format_branch_label = branch_labels 
 
    # options for displaying label colors. 
    if label_colors: 
        if callable(label_colors): 
            def get_label_color(label): 
                return label_colors(label) 
        else: 
            # label_colors is presumed to be a dict 
            def get_label_color(label): 
                return label_colors.get(label, 'black') 
    else: 
        def get_label_color(label): 
            # if label_colors is not specified, use black 
            return 'black' 
 
    # Layout 
 
    def get_x_positions(tree, cladogram=True): 
        """Create a mapping of each clade to its horizontal position. 
 
        Dict of {clade: x-coord} 
        """ 
        global treeDepth
        if cladogram:
            treeDepth = 0
            for terminal in tree.get_terminals():
                candidate = len(tree.get_path(terminal))
                if candidate > treeDepth:
                    treeDepth = candidate
            depths = tree.depths(unit_branch_lengths=True)
            for clade in depths.keys():
                if clade.is_terminal():
                    depths[clade] = treeDepth
        else:
            depths = tree.depths() 
            # If there are no branch lengths, assume unit branch lengths 
            if not max(depths.values()): 
                depths = tree.depths(unit_branch_lengths=True) 
        return depths 
 
    def get_y_positions(tree): 
        global heights
        """Create a mapping of each clade to its vertical position. 
 
        Dict of {clade: y-coord}. 
        Coordinates are negative, and integers for tips. 
        """ 
        maxheight = 2*np.pi
        # Rows are defined by the tips 
        heights = {tip: maxheight*(1 - (i/tree.count_terminals())) 
                   for i, tip in enumerate(reversed(tree.get_terminals()))}
 
        # Internal nodes: place at midpoint of children 
        def calc_row(clade): 
            for subclade in clade: 
                if subclade not in heights: 
                    calc_row(subclade) 
            # Closure over heights 
            heights[clade] = (heights[clade.clades[0]] + 
                              heights[clade.clades[-1]]) / 2.0 
 
        if tree.root.clades: 
            calc_row(tree.root) 
        return heights 
 
    x_posns = get_x_positions(tree) 
    y_posns = get_y_positions(tree) 
    # The function draw_clade closes over the axes object 
    if axes is None: 
        fig = plt.figure() 
        axes = fig.add_subplot(111, projection='polar')
    elif not isinstance(axes, plt.matplotlib.axes.Axes): 
        raise ValueError("Invalid argument for axes: %s" % axes) 
 
    def draw_clade_lines(use_linecollection=False, orientation='horizontal', 
                         y_here=0, x_start=0, x_here=0, y_bot=0, y_top=0, 
                         color=lineColor, lw='.1', resolution=100): 
        """Create a line with or without a line collection object. 
 
        Graphical formatting of the lines representing clades in the plot can be 
        customized by altering this function. 
        """ 
        if orientation == 'horizontal': 
            if y_here == np.pi/2:
                axes.plot([y_here+0.00018, y_here+0.00018], 
                          [x_start+0.012, x_here+0.012], 
                          color=color, lw=lw)
            else:
                axes.plot([y_here, y_here], 
                          [x_start, x_here], 
                          color=color, lw=lw)
        elif orientation == 'vertical': 
            axes.plot(np.linspace(y_bot, y_top, resolution), [x_here]*resolution, color=color, lw=lw)
 
    def draw_clade(clade, x_start, color, lw): 
        """Recursively draw a tree, down from the given clade.""" 
        x_here = x_posns[clade] 
        y_here = y_posns[clade] 
        # phyloXML-only graphics annotations 
        if hasattr(clade, 'color') and clade.color is not None: 
            color = clade.color.to_hex() 
        if hasattr(clade, 'width') and clade.width is not None: 
            lw = clade.width * plt.rcParams['lines.linewidth'] 
        # Draw a horizontal line from start to here 
        draw_clade_lines(use_linecollection=False, orientation='horizontal', 
                         y_here=y_here, x_start=x_start, x_here=x_here, color=color, lw=lw) 
        # Add node/taxon labels 
        '''
        label = label_func(clade) 
        if label not in (None, clade.__class__.__name__): 
            if label in mark:
                axes.text(y_here, x_here, ' %s' % 
                          label, verticalalignment='center', 
                          color=markColor, fontweight=markWeight,
                          rotation=y_here + np.pi/2) 
            else:
                axes.text(y_here, x_here, ' %s' % 
                          label, verticalalignment='center', 
                          color=get_label_color(label),
                          rotation=y_here + np.pi/2) 
        '''
        # Add label above the branch (optional) 
        conf_label = format_branch_label(clade) 
        if conf_label:
            axes.text(0.5 * (x_start + x_here), y_here, conf_label, 
                      fontsize='small', horizontalalignment='center') 
        if clade.clades: 
            # Draw a vertical line connecting all children 
            y_top = y_posns[clade.clades[0]] 
            y_bot = y_posns[clade.clades[-1]] 
            # Only apply widths to horizontal lines, like Archaeopteryx 
            draw_clade_lines(use_linecollection=False, orientation='vertical', 
                             x_here=x_here, y_bot=y_bot, y_top=y_top, 
                             color=color, lw=lw) 
            # Draw descendents 
            for child in clade: 
                draw_clade(child, x_here, color, lw) 
 
    draw_clade(tree.root, 0, 'k', plt.rcParams['lines.linewidth']) 
 
    # If line collections were used to create clade lines, here they are added 
    # to the pyplot plot. 
    for i in horizontal_linecollections: 
        axes.add_collection(i) 
    for i in vertical_linecollections: 
        axes.add_collection(i) 
 
    # Aesthetics 
    
    #xmax = max(x_posns.values()) 
    #axes.set_xlim(-0.05 * xmax, 1.25 * xmax) 
    axes.set_rmax(treeDepth + 1)
    '''
    ticks = np.linspace(0, 2*np.pi, 1 + len(tree.get_terminals()))
    axes.set_xticks(ticks)
    axes.set_xticklabels(leaf.name for leaf in tree.get_terminals())

    angles = np.linspace(0,2*np.pi,len(axes.get_xticklabels())+1)
    #angles[np.cos(angles) < 0] = angles[np.cos(angles) < 0] + np.pi
    angles = np.rad2deg(angles)
    
    fig.canvas.draw()
    labels = []
    for label, theta in zip(axes.get_xticklabels(), angles):
        lab = axes.text(treeDepth, theta, label.get_text(), 
                      transform=label.get_transform(),
                      ha=label.get_ha(), va=label.get_va())
        if np.cos(theta) < 0:
            lab.set_rotation(theta + np.pi)
        else: lab.set_rotation(theta)
        labels.append(lab)
    axes.set_xticklabels([])
    '''
    # Parse and process key word arguments as pyplot options 
    for key, value in kwargs.items(): 
        try: 
            # Check that the pyplot option input is iterable, as required 
            [i for i in value] 
        except TypeError: 
            raise ValueError('Keyword argument "%s=%s" is not in the format ' 
                             'pyplot_option_name=(tuple), pyplot_option_name=(tuple, dict),' 
                             ' or pyplot_option_name=(dict) ' 
                             % (key, value)) 
        if isinstance(value, dict): 
            getattr(plt, str(key))(**dict(value)) 
        elif not (isinstance(value[0], tuple)): 
            getattr(plt, str(key))(*value) 
        elif (isinstance(value[0], tuple)): 
            getattr(plt, str(key))(*value[0], **dict(value[1])) 
 
    if do_show: 
        plt.show()

def drawTree(treeName,  
             outX,
             outY=10,
             outFormat='png', 
             outDPI=300,
             ladderize=False, 
             fontSize=5,
             lineWidth=0.25,
             mark=[]):

    plt.rc('font', size=fontSize)
    plt.rc('lines', lw=lineWidth, color='k')
    plt.rc('figure', figsize=(outX,10))

    tree = Phylo.read(treeName + '.nwk', "newick")
    if ladderize: tree.ladderize()
    protein7 = tree.common_ancestor(toMark[0], toMark[1])
    protein7.color = 'r'

    # for m in toMark:
        # if tree.find_any(m):
            # print('qwe')
            # tree.find_any(m).color = 'r'

    drawMark(tree, lambda n: n.name, do_show=False, mark=[])
    pylab.axis("off")
    pylab.savefig("{0}.{1}".format(treeName, outFormat),
                  format=outFormat,
                  #bbox_inches='tight',
                  #pad_inches=0,
                  dpi=outDPI)

toMark = ['XP_021076601.1_H(+)/Cl(-)_exchange_transporter_7_isoform_X3_Mus_pahari',
          'XP_012063553.1_PREDICTED:_H(+)/Cl(-)_exchange_transporter_7_Atta_cephalotes']
drawTree('1000', outX=10, mark=toMark)