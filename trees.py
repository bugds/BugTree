#!/usr/bin/env python3

from Bio import Phylo
import matplotlib as plt
import pylab

def drawMark(tree, label_func=str, do_show=True, show_confidence=True, 
             axes=None, branch_labels=None, label_colors=None, mark=[],
             markColor='r', markWeight='bold', *args, **kwargs):
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
 
    def get_x_positions(tree): 
        """Create a mapping of each clade to its horizontal position. 
 
        Dict of {clade: x-coord} 
        """ 
        depths = tree.depths() 
        # If there are no branch lengths, assume unit branch lengths 
        if not max(depths.values()): 
            depths = tree.depths(unit_branch_lengths=True) 
        return depths 
 
    def get_y_positions(tree): 
        """Create a mapping of each clade to its vertical position. 
 
        Dict of {clade: y-coord}. 
        Coordinates are negative, and integers for tips. 
        """ 
        maxheight = tree.count_terminals() 
        # Rows are defined by the tips 
        heights = {tip: maxheight - i 
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
        axes = fig.add_subplot(1, 1, 1) 
    elif not isinstance(axes, plt.matplotlib.axes.Axes): 
        raise ValueError("Invalid argument for axes: %s" % axes) 
 
    def draw_clade_lines(use_linecollection=False, orientation='horizontal', 
                         y_here=0, x_start=0, x_here=0, y_bot=0, y_top=0, 
                         color='black', lw='.1'): 
        """Create a line with or without a line collection object. 
 
        Graphical formatting of the lines representing clades in the plot can be 
        customized by altering this function. 
        """ 
        if not use_linecollection and orientation == 'horizontal': 
            axes.hlines(y_here, x_start, x_here, color=color, lw=lw) 
        elif use_linecollection and orientation == 'horizontal': 
            horizontal_linecollections.append(mpcollections.LineCollection( 
                [[(x_start, y_here), (x_here, y_here)]], color=color, lw=lw),) 
        elif not use_linecollection and orientation == 'vertical': 
            axes.vlines(x_here, y_bot, y_top, color=color) 
        elif use_linecollection and orientation == 'vertical': 
            vertical_linecollections.append(mpcollections.LineCollection( 
                [[(x_here, y_bot), (x_here, y_top)]], color=color, lw=lw),) 
 
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
        draw_clade_lines(use_linecollection=True, orientation='horizontal', 
                         y_here=y_here, x_start=x_start, x_here=x_here, color=color, lw=lw) 
        # Add node/taxon labels 
        label = label_func(clade) 
        if label not in (None, clade.__class__.__name__): 
            if label in mark:
                axes.text(x_here, y_here, ' %s' % 
                          label, verticalalignment='center', 
                          color=markColor, fontweight=markWeight) 
            else:
                axes.text(x_here, y_here, ' %s' % 
                          label, verticalalignment='center', 
                          color=get_label_color(label)) 
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
            draw_clade_lines(use_linecollection=True, orientation='vertical', 
                             x_here=x_here, y_bot=y_bot, y_top=y_top, color=color, lw=lw) 
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
 
    if hasattr(tree, 'name') and tree.name: 
        axes.set_title(tree.name) 
    axes.set_xlabel('branch length') 
    axes.set_ylabel('taxa') 
    # Add margins around the tree to prevent overlapping the axes 
    xmax = max(x_posns.values()) 
    axes.set_xlim(-0.05 * xmax, 1.25 * xmax) 
    # Also invert the y-axis (origin at the top) 
    # Add a small vertical margin, but avoid including 0 and N+1 on the y axis 
    axes.set_ylim(max(y_posns.values()) + 0.8, 0.2) 
 
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
             outDPI=100,
             ladderize=True, 
             fontSize=12,
             lineWidth=0.75,
             mark=[]):

    # ax = plt.pyplot.subplot(111, projection='polar')

    plt.rc('font', size=fontSize)
    plt.rc('lines', lw=lineWidth, color='k')
    plt.rc('figure', figsize=(outX,10))

    tree = Phylo.read(treeName + '.nwk', "newick")
    if ladderize: tree.ladderize()

    # for m in toMark:
        # if tree.find_any(m):
            # print('qwe')
            # tree.find_any(m).color = 'r'

    drawMark(tree, lambda n: n.name, do_show=False, mark=mark)
    pylab.axis("off")
    pylab.savefig("{0}.{1}".format(treeName, outFormat),
                  format=outFormat,
                  bbox_inches='tight',
                  pad_inches=0,
                  dpi=outDPI)

toMark = ['NP_001278.1 Homo sapiens (human)',
          'NP_001278.1 H/Cl-exchange transporter 7']
drawTree('orthologs', outX=12, mark=toMark, outFormat='svg')
drawTree('orthologs', outX=10, mark=toMark)
drawTree('paralogs', outX=12, mark=toMark, outFormat='svg')
drawTree('paralogs', outX=10, mark=toMark)

# TODO:
# 1) Circular plot