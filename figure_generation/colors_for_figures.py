
# A color blind/friendly color cycle for Matplotlib line plots
# from https://gist.github.com/thriveth/8560036
color_blind_friendly_color_cycle_analogs = {'blue': '#377eb8', 'orange': '#ff7f00', 'green': '#4daf4a',
                                            'pink': '#f781bf', 'brown': '#a65628', 'purple': '#984ea3',
                                            'gray': '#999999', 'red': '#e41a1c', 'yellow': '#dede00'}




#from  https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib
#from, https://gist.github.com/ihincks/6a420b599f43fcd7dbd79d56798c4e5a
def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    result= colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])
    positive_r=[max([0,r])for r in result]
    print(result)
    return positive_r