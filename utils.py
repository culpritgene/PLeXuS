from imports import *

def draw_figure(canvas, figure, loc=(0, 0)):
    """ Draw a matplotlib figure onto a Tk canvas

    loc: location of top-left corner of figure on canvas in pixels.
    Inspired by matplotlib source: lib/matplotlib/backends/backend_tkagg.py

    """
    figure_canvas_agg = FigureCanvasAgg(figure)
    figure_canvas_agg.draw()
    figure_x, figure_y, figure_w, figure_h = figure.bbox.bounds
    figure_w, figure_h = int(figure_w), int(figure_h)
    photo = tk.PhotoImage(master=canvas, width=figure_w, height=figure_h)

    # Position: convert from top-left anchor to center anchor
    canvas.create_image(loc[0] + figure_w/2, loc[1] + figure_h/2, image=photo)

    # Unfortunately, there's no accessor for the pointer to the native renderer
    tkagg.blit(photo, figure_canvas_agg.get_renderer()._renderer, colormode=2)

    # Return a handle which contains a reference to the photo object
    # which must be kept live or else the picture disappears
    return photo

class CONFIG():
    expl = 'Just a class to hold config values and save them with the network'
    def __init__(self):
        pass

def get_uniprot_data(uniprot_name):
    handle = requests.get(f"http://www.uniprot.org/uniprot/{uniprot_name}.xml", headers=header)
    with open('TMP', 'w') as tmp:
        tmp.write(handle.text)
    record = SeqIO.read('TMP', "uniprot-xml")
    #os.remove('TMP')
    return record

def print_func_annotations(*args, func='comment_function'):
    for v in args:
        try:
            print(v.annotations[func])
            print('#############################'*3)
        except IndexError:
            pass

header = {'User-Agent':'python culpritgene@gmail.com'}