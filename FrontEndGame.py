from imports import *
from utils import *
from BackEndGame import *

class Gui():
    def __init__(self, root, init_pa):
        self.root = root
        self.init_pa = init_pa
        self.prot_names = {}
        self.prot_names.update({'kegg': js.load(open('Resources/Kegg_prot_Names.json', 'r')),
                                'signor': js.load(open('Resources/Signor_pathways.json', 'r')),
                                'panther': js.load(open('Resources/Panther_db_selected.json', 'r')),
                                'reactome': []})

        self.init_network_filt = self.init_pa.graph
        self.current_game_network = None
        self.current_labirinth_game = None
        self.current_game_network_filt = None
        self.known_locations = {}

        self.graph1 = False

        self.CANVAS_WIDTH = 800
        self.CANVAS_HEIGHT = 1100
        self.entry = tk.Entry(root)
        self.USER = tk.StringVar()
        self.USER.set("Unknown")
        self.GMODE = tk.StringVar()
        self.GMODE.set("Labirinth Game")
        self.LMODE = tk.StringVar()
        self.LMODE.set("Autogenerate")

        self.CURSOR_POSITION = tk.StringVar()
        self.CURSOR_POSITION.set("X:0 | Y:0")
        self.CURR_INFO = tk.StringVar()
        self.CURR_INFO.set('')
        self.GAME_ON = False

        self.DB_LIST = ['kegg', 'signor', 'panther', 'reactome']
        self.CUR_DB = tk.StringVar()
        self.CUR_DB.set("kegg")

        self.SELECTED_GENE_SET = tk.StringVar()

        self.PATHWAY_LIST = list(self.prot_names[self.CUR_DB.get()])



        self.CUR_NODES = ['random']
        self.CUR_F_NODE = tk.StringVar()
        self.CUR_F_NODE.set('random')

        self.CUR_PATHWAY = tk.StringVar()
        self.CUR_PATHWAY.set("None")
        self.CURR_SUBNET_SIZE = tk.StringVar()
        self.CURR_SUBNET_SIZE.set('Na')
        self.CURR_ININET_SIZE = tk.StringVar()
        self.CURR_ININET_SIZE.set('10080')
        self.CURR_STEP = tk.IntVar()
        self.CURR_STEP.set(0)

        self.FILTER_REF_NUM_sele = 0
        self.FILTER_DEGREE_sele = 0
        self.FILTER_REF_NUM_init = 0
        self.FILTER_DEGREE_init = 0

        self.LAYOUT_type = tk.StringVar()
        self.LAYOUT_type.set('reing')
        self.FILTER_TYPE = tk.StringVar()
        self.FILTER_TYPE.set('Ref_number')

        self.REMOVE_SELF_CONN = tk.BooleanVar()
        self.REMOVE_SELF_CONN.set(False)
        self.SELECTION_TYPE = tk.StringVar()
        self.SELECTION_TYPE.set('random walk')
        self.MASK_VOL = tk.StringVar()
        self.MASK_VOL.set('subgraph')
        self.TOTAL_LIVES = tk.IntVar()
        self.TOTAL_LIVES.set(2)
        self.CONTRASTIVE_STEP = 2

        self.LABIRINTH_CONFIG = CONFIG()
        self.LABIRINTH_CONFIG.STATS_COEFF = tk.DoubleVar()
        self.LABIRINTH_CONFIG.STATS_COEFF.set(0.1)
        self.LABIRINTH_CONFIG.MAX_LENGTH = tk.IntVar()
        self.LABIRINTH_CONFIG.MAX_LENGTH.set(2)
        self.LABIRINTH_CONFIG.BRANCHING_NUM = tk.IntVar()
        self.LABIRINTH_CONFIG.BRANCHING_NUM.set(3)
        self.LABIRINTH_CONFIG.CONNECTGRAPH = BooleanVar()
        self.LABIRINTH_CONFIG.CONNECTGRAPH.set(False)
        self.LABIRINTH_CONFIG.GAME_MODE = 'prime'

        #    self.LABIRINTH_CONFIG.FILTER_DEGREE = 0

        self.STATS = CONFIG()
        self.STATS.PLAYED = tk.IntVar(0)
        self.STATS.WINLOS = tk.StringVar()
        self.STATS.WINLOS.set('0/0')
        self.STATS.KNOWN = tk.IntVar(0)
        self.STATS.stats_proteins = np.array([0])
        self.STATS.stats_games = np.array([0,0,0])

        self.known_locations = {}
        self.TEXT_HIGHLIGHT = tk.StringVar()
        self.TEXT_HIGHLIGHT.set('')

        self.canvas = tk.Canvas(self.root, width=self.CANVAS_WIDTH, height=self.CANVAS_HEIGHT, background='white')
        self.canvas.grid(row=1, column=2)

        self.turtle = turtle.RawTurtle(self.canvas)
        self.turtle.pensize(width=2)
        #    self.turtle.hideturtle()

        frame0 = Frame(self.root)  ## Here we create a separate frame1 within the class Gui!
        frame0.grid(row=0, column=2, sticky="n")
        self.frame0 = frame0

        frame01 = Frame(self.root, width=250, height=1100)  ## Here we create a separate frame1 within the class Gui!
        frame01.grid(row=1, column=0, sticky="n", )
        frame01.grid_propagate(False)
        self.frame01 = frame01

        frame1 = Frame(frame01, width=250, height=350)  ## Here we create a separate frame1 within the class Gui!
        frame1.grid(row=1, column=0, sticky="n", )
        frame1.grid_propagate(False)
        self.frame1 = frame1

        frame1opt = Frame(frame01, width=250, height=350)  ## Frame to hold options parameters!
        frame1opt.grid(row=1, column=0, sticky="n", )
        frame1opt.grid_propagate(False)
        self.frame1opt = frame1opt

        frame1stats = Frame(frame01, width=250, height=350)  ## Frame to hold options parameters!
        frame1stats.grid(row=1, column=0, sticky="n", )
        frame1stats.grid_propagate(False)
        self.frame1stats = frame1stats

        frame2 = Frame(frame01, width=250, height=750)  ## Here we create a separate frame1 within the class Gui!
        frame2.grid(row=2, column=0, sticky="S", )
        frame2.grid_propagate(False)

        self.option = tk.OptionMenu(frame1, self.USER, "Culpritgene", "Unknown", "Other",
                                    command=self.change_user_command)
        label1 = Label(frame1, text="User").grid(row=0, column=0, sticky="nw", pady=7)
        label2 = Label(frame1, text="Game Mode").grid(row=1, column=0, sticky="w")
        label3 = Label(frame1, text="Learning Mode").grid(row=2, column=0, sticky="w")
        self.option.grid(row=0, column=1, sticky="nwe")
        self.game_mode = tk.OptionMenu(frame1, self.GMODE, "Labirinth Game", 'Blocks Game')
        self.learning_mode = tk.OptionMenu(frame1, self.LMODE, "Autogenerate", "User Pick", "Load Tasks")
        self.game_mode.grid(row=1, column=1, sticky="nwe")
        self.learning_mode.grid(row=2, column=1, sticky="nwe")
        # entry = Entry(frame1).grid(row = 1,column = 1,sticky = E+ W)
        # entry1 = Entry(frame1).grid(row = 2,column = 1, sticky = E)
        self.Options_B = Button(frame1, command=self.options_command, text="Options").grid(row=3, column=1, sticky="we")
        self.Stats_B = Button(frame1, text="Stats", command=self.stats__command).grid(row=4, column=1, sticky="we")
        self.About_B = Button(frame1, text="About").grid(row=5, column=1, sticky="we")
        self.Draw_B = Button(frame1, command=self.Show_uniprot2, text="Get New Hint").grid(row=6, column=1, sticky="we")
        self.Show_B = Button(frame1, command=self.show_WHOLE_graph_command, text="Draw Graph").grid(row=7, column=1,
                                                                                                    sticky="we")
        self.Receite_B = Button(frame1, text="Receite",  command=self.receite).grid(row=8, column=1, sticky="we",)
        self.Clear_B = Button(frame1, text="Clear", command=self.clear).grid(row=9, column=1, sticky="we")
        self.SAVE = Button(frame1, command=self.save_command, text="Save Stats").grid(row=10, column=1, sticky="we")
        self.Exit_B = Button(frame1, command=self.exit_command, text="Exit").grid(row=11, column=1, sticky="we")

        ################################################ OPTIONS
        label1opt = Label(frame1opt, text="Game Options").grid(row=0, column=1, sticky="nw", pady=8)

        GameOptFrame = Frame(frame1opt)
        GameOptFrame.grid(row=1, column=1, sticky="nw", pady=5)
        SelectionOptFrame = Frame(frame1opt)
        SelectionOptFrame.grid(row=2, column=1, sticky="nw", pady=5)
        HintsOptFrame = Frame(frame1opt)
        HintsOptFrame.grid(row=3, column=1, sticky="nw", pady=5)
        DisplayOptFrame = Frame(frame1opt)
        DisplayOptFrame.grid(row=4, column=1, sticky="nw", pady=5)

        label21opt = Label(GameOptFrame, text="Lives")
        label21opt.grid(row=1, column=0, sticky="w")
        label21opt.config(font=("Arial", 10))
        lives = Entry(GameOptFrame, width=6, text=self.TOTAL_LIVES, textvariable=self.TOTAL_LIVES)
        lives.grid(row=1, column=1, sticky="n", padx=10)
        lives.bind('<Return>', self.set_lives_command)

        label3opt = Label(GameOptFrame, text="Branch")
        label3opt.grid(row=2, column=0, sticky="w")
        label3opt.config(font=("Arial", 10))
        branching = Entry(GameOptFrame, width=6, text=self.LABIRINTH_CONFIG.BRANCHING_NUM,
                          textvariable=self.LABIRINTH_CONFIG.BRANCHING_NUM, )
        branching.grid(row=2, column=1, sticky="n", padx=6)
        branching.bind('<Return>', self.set_branchin_num_command)

        label3opt = Label(GameOptFrame, text="Length")
        label3opt.grid(row=2, column=3, sticky="w")
        label3opt.config(font=("Arial", 10))
        length = Entry(GameOptFrame, width=6, text=self.LABIRINTH_CONFIG.MAX_LENGTH,
                       textvariable=self.LABIRINTH_CONFIG.MAX_LENGTH, )
        length.grid(row=2, column=4, sticky="n", padx=6)
        length.bind('<Return>', self.set_maxlen_command)

        label4opt = Label(SelectionOptFrame, text="Selecting Type")
        label4opt.grid(row=2, column=0, sticky="w")
        label4opt.config(font=("Arial", 10))
        self.SelectionType = tk.OptionMenu(SelectionOptFrame, self.SELECTION_TYPE, "random walk", 'clusters', command=self.change_sele_type_command)
        self.SelectionType.grid(row=2, column=1, sticky="s", padx=4)

        label4opt = Label(SelectionOptFrame, text="Contrastive Step")
        label4opt.grid(row=3, column=0, sticky="w")
        label4opt.config(font=("Arial", 10))
        self.ContrastValue = tk.OptionMenu(SelectionOptFrame, self.CONTRASTIVE_STEP, "unconnected",
                                           "max distance", 'med distance', "closest wrong", )
        self.ContrastValue.grid(row=3, column=1, sticky="s", padx=4)

        label5opt = Label(SelectionOptFrame, text="Filtering Type")
        label5opt.grid(row=4, column=0, sticky="w")
        label5opt.config(font=("Arial", 10))
        self.FilteringType = tk.OptionMenu(SelectionOptFrame, self.FILTER_TYPE, "Ref_number", "degree",
                                           'betwenness_z', 'combined', 'user_score',  command=self.change_FilteringType_command)
        self.FilteringType.grid(row=4, column=1, sticky="s", padx=4)

        label6opt = Label(HintsOptFrame, text="Text Mask degree")
        label6opt.grid(row=2, column=0, sticky="w")
        label6opt.config(font=("Arial", 10))
        self.MaskDegree = tk.OptionMenu(HintsOptFrame, self.MASK_VOL, "subgraph", "neighbors", 'all', command=self.change_masking_command)
        self.MaskDegree.grid(row=2, column=1, sticky="s", padx=4)

        label7opt = Label(DisplayOptFrame, text="Layout Type")
        label7opt.grid(row=2, column=0, sticky="w")
        label7opt.config(font=("Arial", 10))
        self.ContrastValue = tk.OptionMenu(DisplayOptFrame, self.LAYOUT_type, "reing", "trees", "sugi", 'circle', 'kawaii', command=self.change_layout_command)
        self.ContrastValue.grid(row=2, column=1, sticky="w", padx=4)

        label7opt = Label(DisplayOptFrame, text="Connected Graph")
        label7opt.grid(row=3, column=0, sticky="w")
        label7opt.config(font=("Arial", 10))
        self.ConnectGraph = Checkbutton(DisplayOptFrame, variable=self.LABIRINTH_CONFIG.CONNECTGRAPH, onvalue=True, offvalue=False)
        self.ConnectGraph.grid(row=3, column=1, sticky="w", padx=4)

        label22opt = Label(DisplayOptFrame, text="Remove self-edges")
        label22opt.grid(row=4, column=0, sticky="w")
        label22opt.config(font=("Arial", 10))
        lives = Checkbutton(DisplayOptFrame, variable=self.REMOVE_SELF_CONN, onvalue=True, offvalue=False)
        lives.grid(row=4, column=1, sticky="w", padx=4)

        self.OPTIONS_EXIT = Button(frame1opt, command=self.options_exit, text='Return').grid(row=5, column=1,sticky="w", padx=4)
        ##############################################################

        ################################################## STATS
        label1stats = Label(frame1stats, text=f"STATS FOR USER: {self.USER.get()}",
                            textvariable=self.USER).grid(row=0, column=0, sticky="nw", pady=10)

        label2stats = Label(frame1stats, text="Games Played")
        label2stats.grid(row=1, column=0, sticky="w")
        label2stats.config(font=("Arial", 10))
        label22stats = Label(frame1stats, text=self.STATS.PLAYED, textvariable=self.STATS.PLAYED)
        label22stats.grid(row=1, column=1, sticky="w", padx=4)
        label22stats.config(font=("Arial", 10))

        label3stats = Label(frame1stats, text="Win/Loss")
        label3stats.grid(row=2, column=0, sticky="w")
        label3stats.config(font=("Arial", 10))
        label32stats = Label(frame1stats, text=self.STATS.WINLOS, textvariable=self.STATS.WINLOS)
        label32stats.grid(row=2, column=1, sticky="w", padx=4)
        label32stats.config(font=("Arial", 10))

        label4stats = Label(frame1stats, text="Known Proteins")
        label4stats.grid(row=3, column=0, sticky="w")
        label4stats.config(font=("Arial", 10))
        label42stats = Label(frame1stats, text=self.STATS.KNOWN, textvariable=self.STATS.KNOWN)
        label42stats.grid(row=3, column=1, sticky="w", padx=4)
        label42stats.config(font=("Arial", 10))

        label5stats = Label(frame1stats, text="Display tails")
        label5stats.grid(row=4, column=0, sticky="w")
        label5stats.config(font=("Arial", 10))
        button5stats = Button(frame1stats, text='Show Best and Worst', command=self.make_best_and_worst, width=15)
        button5stats.grid(row=4, column=1, sticky="w", padx=4)
        button5stats.config(font=("Arial", 10))

        self.stats_EXIT = Button(frame1stats, command=self.options_exit, text='Return').grid(row=5, column=1, sticky="w", padx=4)
        #############################################################
        # len(np.where(self.get_stats()[:,-1]!=0)[0]))

        ### raise main menu frame back on top
        self.frame1.tkraise()

        ### Buttons for Pathway Selection

        self.Load_Gene_Set = Entry(frame0, width=10, text='', textvariable=self.SELECTED_GENE_SET)
        self.Load_Gene_Set.grid(row=0, column=1, sticky='WN', pady=4, padx=6)
        self.Load_Gene_Set.bind('<Return>', self.load_gene_set_function)

        self.Select_DB = tk.OptionMenu(frame0, self.CUR_DB, *self.DB_LIST, command=self.change_db_command, )
        self.Select_DB.grid(row=0, column=2, sticky="WN")

        self.Select_Pathway_M = tk.OptionMenu(frame0, self.CUR_PATHWAY, *self.PATHWAY_LIST, command=self.Load_network_from_pathway)
        self.Select_Pathway_M.grid(row=0, column=3, sticky="N")

        self.Select_First_Node = tk.OptionMenu(frame0, self.CUR_F_NODE, *self.CUR_NODES, command=self.select_first_node)
        self.Select_First_Node.grid(row=0, column=4, sticky="N")

        self.START_B = Button(frame0, command=self.Start_command, text='Start!').grid(row=0, column=5, sticky="EN")

        self.FILT_ININET_SIZE = Label(frame0, text=self.FILTER_TYPE, textvariable=self.FILTER_TYPE, background='white').grid(row=1, column=1,sticky="SW",pady=3.4, padx=10)

        self.filter_ref_slider_sele = Scale(frame0, orient='horizontal', from_=0, to=200, resolution=10, command=self.set_Ref_filter_value_sele)
        self.filter_ref_slider_sele.grid(row=1,column=2,sticky="S")
        self.filter_ref_slider_init = Scale(frame0, orient='horizontal', from_=0, to=200, resolution=10, command=self.set_Ref_filter_value_init)
        self.filter_ref_slider_init.grid(row=1,column=3,sticky="S")

        self.FILT_NET_SIZE = Label(frame0, textvariable=self.CURR_SUBNET_SIZE,
                                   text=self.CURR_SUBNET_SIZE, background='white').grid(row=1, column=4, sticky="SE",
                                                                                        pady=3.4, padx=5)
        self.FILT_ININET_SIZE = Label(frame0, textvariable=self.CURR_ININET_SIZE, text=self.CURR_ININET_SIZE,
                                      background='white').grid(row=1, column=5, sticky="SE", pady=3.4, padx=5)

        self.TEXT_BOX = Label(frame2, textvariable=self.TEXT_HIGHLIGHT, text=self.TEXT_HIGHLIGHT,
                              background='white').grid(row=10, column=0, sticky="WS", pady=60, padx=15)

        cursor = Label(self.canvas, textvariable=self.CURSOR_POSITION, text=self.CURSOR_POSITION, background='white')
        self.cursor_position = self.canvas.create_window(275, 475, window=cursor)

        self.info = Label(self.canvas, textvariable=self.CURR_INFO, text=self.CURR_INFO, background='white')
        self.info.configure(font=("Arial", 15))
        self.info_box = self.canvas.create_window(-50, 475, window=self.info)

        graph_saver = Button(self.canvas, text='Save Game', command=self.save_game_graph, width=6)
        self.graph_saver = self.canvas.create_window(280, -500, window=graph_saver)
        graph_loader = Button(self.canvas, text='Load Game', command=self.load_game_graph, width=6)
        self.graph_saver = self.canvas.create_window(350, -500, window=graph_loader)

        self.canvas.bind('<Motion>', self.motion)
        self.canvas.bind("<Button-1>", self.user_decision)

    def receite(self):
        msg = 'One round of receiting!'
        if self.USER.get()=='Unknown':
            self.INFO_update('Only registered user can receit!')
        else:
            self.INFO_update(msg, 1200)
            df = self.STATS.stats_proteins
            df2 = df[df['Score']!=0]
            genes = np.random.choice(df2.index, size=20, replace=False)
            pc_proteins = [v.index for v in self.init_pa.graph.vs() if v['label'] in genes]
            pc_subgraph = self.init_pa.graph.induced_subgraph(pc_proteins)
            self.current_game_network = pc_subgraph
            self.current_game_network_filt = pc_subgraph
            self.CURR_SUBNET_SIZE.set(str(len(pc_proteins)))
            self.set_curr_nodes_list()
            self.Start_command()


    def get_stats(self):
        try:
            return np.stack(self.init_pa.graph.vs()[self.USER.get()])
        except:
            return np.zeros(shape=[10, 3])

    def motion(self, event):
        x, y = event.x, event.y
        self.CURSOR_POSITION.set('X:' + str(x) + '|' + 'Y:' + str(y))
        r = 18
        enter = 0
        for lab in list(self.known_locations):
            c, tick = self.known_locations[lab]
            x0, y0, x1, y1 = c.flatten()
            if self.adjx(x1) - r < x < self.adjx(x1) + r and self.adjy(y1) - r < y < self.adjy(y1) + r:
                enter += 1
                annot = self.current_labirinth_game.PROT_INFO[lab]
                annot = self.split_by_ns(annot, add=lab)
                self.TEXT_HIGHLIGHT.set(annot)
        if enter == 0:
            self.TEXT_HIGHLIGHT.set('')
        if enter > 1:
            self.CURR_INFO.set('<probable nodes overlapp...>')
        else:
            pass
        # self.CURR_INFO.set('')

    def user_decision(self, event):
        if self.GAME_ON:
            x, y = event.x, event.y
            r = 18
            zones = []
            for lab in list(self.known_locations):
                c, tick = self.known_locations[lab]
                x0, y0, x1, y1 = c.flatten()
                if tick == 0:
                    if self.adjx(x1) - r < x < self.adjx(x1) + r and self.adjy(y1) - r < y < self.adjy(y1) + r:
                        corr = self.current_labirinth_game.true_path[self.CURR_STEP.get() + 1]
                        self.current_labirinth_game.path_chosen.append([lab, corr])
                        if lab == corr:
                            #### CORRECT DECISION
                            msg = 'Correct! Next step!'
                            self.CURR_INFO.set(msg)
                            self.turtle_run(c, lab, line='green', text='green', speed=5)
                            for l in self.current_labirinth_game.game_subgraph.vs()['label'][1:][
                                     self.LABIRINTH_CONFIG.BRANCHING_NUM.get() * self.CURR_STEP.get():self.LABIRINTH_CONFIG.BRANCHING_NUM.get() * (self.CURR_STEP.get() + 1)]:
                                self.known_locations[l][-1] += 1
                                if l != lab:
                                    c, tickk = self.known_locations[l]
                                    self.turtle_run(c, l, line='red', text='red', speed=7)
                            print(msg)
                            self.CURR_STEP.set(self.CURR_STEP.get() + 1)
                            if lab == self.current_labirinth_game.true_path[-1] and self.CURR_STEP.get() == len(
                                    self.current_labirinth_game.true_path) - 1:
                                msg = f'HOOREY, THIS GAME IS WON!'
                                self.CURR_INFO.set(msg)
                                print(msg)
                                self.GAME_ON = False
                                self.register_stats(res=1)
                            else:
                                self.draw_graph_command()
                        else:
                            #### INCORRECT DECISION
                            self.current_labirinth_game.lives -= 1
                            self.turtle_run(c, lab, line='red', text='red', speed=5)
                            if self.current_labirinth_game.lives > 0:
                                msg = f'Sorry, this time wrong. Lets try again! (lives left:{self.current_labirinth_game.lives})'
                                self.CURR_INFO.set(msg)
                            else:
                                msg = f'Sorry, game over.'
                                self.GAME_ON = False
                                self.CURR_INFO.set(msg)
                                self.register_stats(res=2)


    def show_WHOLE_graph_command(self):

        layout = self.current_game_network_filt.layout_fruchterman_reingold(
            repulserad=self.current_game_network_filt.vcount() ** 2.2,
            maxiter=1000, area=self.current_game_network_filt.vcount() ** 2, coolexp=2, )

        bbox = igraph.BoundingBox(0, 0, self.CANVAS_WIDTH / 2, self.CANVAS_HEIGHT / 2)
        layout.fit_into(bbox)
        layout.center(0, 0)
        edges = []
        for e in self.current_game_network_filt.es():
            edges.append(e.tuple)
        edges = np.array(edges)
        coords = np.array(layout.coords)

        self.crs = list(coords[edges].flatten())
        self.graph1 = self.canvas.create_line(*self.crs, fill="blue", width=2)

    def draw_graph_command(self, type_='new'):
        edges = self.current_labirinth_game.edges
        coords = self.current_labirinth_game.coords
        self.zz = list(zip(coords[edges], self.current_labirinth_game.game_subgraph.vs()['label'][1:]))
        if self.CURR_STEP.get() == 0:
            self.turtle.penup()
            self.turtle.setposition(self.zz[0][0][0])
            self.turtle.pendown()
            self.turtle.pencolor('blue')
            self.turtle.write(self.current_labirinth_game.start_node, font=('Verdana', 16, 'bold'), align='center')
            self.turtle.pencolor('black')
            self.known_locations.update({self.current_labirinth_game.start_node: [self.zz[0][0][[1, 0]], 1]})
        zz_slice = self.zz[self.LABIRINTH_CONFIG.BRANCHING_NUM.get() * self.CURR_STEP.get():self.LABIRINTH_CONFIG.BRANCHING_NUM.get() * (self.CURR_STEP.get() + 1)]
        for c, lab in shuffle(zz_slice):
            self.turtle_run(c, lab)

    def turtle_run(self, c, lab, line='black', text='blue', speed=3):
        self.turtle.pencolor(line)
        self.turtle.speed(speed)
        self.turtle.penup()
        self.turtle.setposition(c[0])
        self.turtle.pendown()
        self.turtle.setposition(c[1])
        if lab in self.known_locations.keys():
            self.known_locations[lab][-1] += 1
        else:
            self.known_locations.update({lab: [c, 0]})
            self.turtle.circle(3)
            self.turtle.pencolor(text)
            self.turtle.write(lab, font=('Verdana', 16, 'bold'), align='center')
            self.turtle.pencolor(line)

    def split_by_ns(self, text, add=None, winsize=30):
        toktext = list(text[2:-2])
        max_iter = len(toktext) / winsize
        for k in range(1, int(max_iter + 1)):
            toktext.insert(winsize * k, '\n')
        toktext.insert(0, '{:{align}{width}}\n'.format('Functional Annotation', align='^', width=str(winsize + 2)))
        toktext.insert(1, '{:{align}{width}}\n'.format(add, align='^', width=str(winsize + 2)))
        toktext = ''.join(toktext)
        ll = toktext.split('\n')[-1]
        if len(ll) > winsize:
            toktext = toktext[:-(len(ll) - winsize)][:-3]
            toktext += '...'
        return toktext


    def INFO_update(self, msg, time=1200):
        prev_msg = self.CURR_INFO.get()
        self.CURR_INFO.set(msg)
        self.info.update()
        self.info.after(time)
        self.CURR_INFO.set(prev_msg)

    def save_game_graph(self):
        graphs = [0]
        for root, dirs, files in os.walk("."):
            for filename in files:
                if filename.split('.')[-1] == 'pa':
                    graphs.append(int(filename.split('game_')[-1][:-3]))
        iterr = max(graphs) + 1
        net = copy.copy(self.current_labirinth_game)
        del net.parent
        pickle.dump(net, open(f'Resources/saved_sessions/game_{iterr}.pa', 'wb'))

    def save_game_graph2(self):
        a,b,c = self.current_labirinth_game.game_subgraph, self.current_labirinth_game.true_path, self.current_labirinth_game.path_variants


    def load_game_graph(self, path='Resources/saved_sessions/game_5.pa'):
        self.GAME_ON = True
        self.CURR_INFO.set('Loading default game net!')
        net = pickle.load(open(path, 'rb'))
        self.current_labirinth_game = net
        self.current_labirinth_game.parent = self
        self.draw_graph_command()
        self.CURR_INFO.set('Choose right leaf!')

    def save_command(self):
        # pickle.dump(self.init_pa, open('Network_Last_Burn.pa','wb'))
        st = copy.copy(self.STATS)
        del st.PLAYED
        del st.WINLOS
        del st.KNOWN
        pickle.dump(st, open(f'Resources/Stats_{self.USER.get()}.st', 'wb'))

    def exit_command(self, ):
        self.root.destroy()

    def change_user_command(self, x):
        self.USER.set(x)
        if self.USER.get() != 'Unknown':
            x = pickle.load(open(f'Resources/Stats_{self.USER.get()}.st', 'rb'))
            self.STATS.stats_games = x.stats_games
            self.STATS.stats_proteins = x.stats_proteins
            self.update_STATS_vars()

    def change_layout_command(self, x):
        self.LAYOUT_type.set(x)

    def update_STATS_vars(self):
        pl, w, l = self.STATS.stats_games
        df = self.STATS.stats_proteins
        self.STATS.PLAYED.set(pl)
        self.STATS.WINLOS.set(f'{w}/{l}')
        self.STATS.KNOWN.set(df[df.sum(axis=1) != 0].shape[0])

    def set_branchin_num_command(self, x):

        self.LABIRINTH_CONFIG.BRANCHING_NUM.set(x.__dict__['widget'].get())

    def set_maxlen_command(self, x):
        self.LABIRINTH_CONFIG.MAX_LENGTH.set(x.__dict__['widget'].get())

    def set_lives_command(self, x):
        self.TOTAL_LIVES.set(x.__dict__['widget'].get())

    def about_command(self, ):
        pass

    def stats__command(self, ):
        self.frame1stats.tkraise()

    def options_command(self, ):
        self.frame1opt.tkraise()

    def options_exit(self):
        self.frame1.tkraise()

    def change_FilteringType_command(self, x):
        self.FILTER_TYPE.set(x)
        m = np.array(self.init_pa.graph.vs()[x])
        min_, max_, resol = m.min(), np.quantile(m, 0.99), (np.quantile(m, 0.99)-m.min())/30
        if x=='Ref_number':
            resol=10
            max_=200
        self.filter_ref_slider_init.configure(from_=min_, to= max_, resolution=resol)
        self.filter_ref_slider_sele.configure(from_=min_, to= max_, resolution=resol)


    def change_sele_type_command(self, x):
        self.SELECTION_TYPE.set(x)

    def change_masking_command(self, x):
        self.MASK_VOL.set(x)

    def register_stats(self, res=1):
        self.STATS.stats_games[0] += 1
        self.STATS.stats_games[res] += 1
        disc = self.STATS.stats_proteins
        labs = disc.index
        assert set(labs) == set(self.init_pa.graph.vs()['label']), print('Mismatch in Protein numbers!')
        for x in self.current_labirinth_game.path_chosen:
            if x[0] == x[1]:
                c = (1 + self.LABIRINTH_CONFIG.STATS_COEFF.get()) ** self.LABIRINTH_CONFIG.BRANCHING_NUM.get()
                disc.at[x[1], 'Pos'] += 1
                disc.at[x[1], 'Score'] += c
            else:
                c = -(1 - self.LABIRINTH_CONFIG.STATS_COEFF.get()) ** self.LABIRINTH_CONFIG.BRANCHING_NUM.get()
                disc.at[x[1], 'Neg'] += 1
                disc.at[x[1], 'Score'] += c
        self.update_STATS_vars()


    def adjx(self, x):
        return x + self.CANVAS_WIDTH / 2

    def adjy(self, y):
        return self.CANVAS_HEIGHT / 2 - y

    def clear(self, ):
        self.turtle.setposition(0, 0)
        self.turtle.clear()
        self.current_labirinth_game = None
        self.CURR_STEP.set(0)
        self.known_locations = {}
        self.TEXT_HIGHLIGHT.set('')
        self.CURR_INFO.set('')
        if self.graph1:
            self.canvas.delete(self.graph1)

    def random_first_node(self):
        start_node_probs = np.array(self.current_game_network_filt.degree())
        start_node_probs[start_node_probs == 0] = 1000000
        start_node_probs = 1 / (start_node_probs ** 2)
        start_node_probs = start_node_probs / (sum(start_node_probs))
        start_node = np.random.choice(self.current_game_network_filt.vs()['label'],  p=start_node_probs)
        return start_node

    def Start_command(self, ):
        self.GAME_ON = True
        start_node = self.CUR_F_NODE.get()
        if start_node=='random':
            start_node = self.random_first_node()
        self.current_labirinth_game = labirinth_game(self, path_length=self.LABIRINTH_CONFIG.MAX_LENGTH.get(),
                                                     mode='prime', start_node=start_node)
        self.current_labirinth_game.start_game_01()
        self.draw_graph_command()
        self.CURR_INFO.set('Choose right leaf!')



    def set_curr_nodes_list(self, x='a'):
        self.CUR_NODES = ['random']+self.current_game_network_filt.vs()['label']
        self.CUR_F_NODE.set('random')
        self.Select_First_Node['menu'].delete(0, 'end')
        for choice in self.CUR_NODES:
            self.Select_First_Node['menu'].add_command(label=choice, command=tk._setit(self.CUR_F_NODE, choice,
                                                                                      self.select_first_node))

    def change_db_command(self, x):
        self.CUR_PATHWAY.set('')
        self.PATHWAY_LIST = list(self.prot_names[x].keys())
        print(x)
        self.Select_Pathway_M['menu'].delete(0, 'end')
        for choice in self.PATHWAY_LIST:
            self.Select_Pathway_M['menu'].add_command(label=choice, command=tk._setit(self.CUR_PATHWAY, choice,
                                                                                      self.Load_network_from_pathway))

    def load_gene_set_function(self, x):
        self.SELECTED_GENE_SET.set(x.__dict__['widget'].get())
        self.root.after(1200)
        genes = self.SELECTED_GENE_SET.get().split()
        self.SELECTED_GENE_SET.set('')
        self.Load_network_from_GS(genes)
        self.set_curr_nodes_list()
        msg = 'Loaded Users Gene Set!'
        self.INFO_update(msg, 2000)

    def select_first_node(self, x):
        self.CUR_F_NODE.set(x)


    def Load_network_from_pathway(self, x):
        self.CUR_PATHWAY.set(x)
        pc_proteins = [i for i in self.prot_names[self.CUR_DB.get()][self.CUR_PATHWAY.get()] if
                       i in self.init_pa.nodDct]
        pc_subgraph = self.init_pa.graph.induced_subgraph(pc_proteins)
        self.current_game_network = pc_subgraph
        self.current_game_network_filt = pc_subgraph
        self.CURR_SUBNET_SIZE.set(str(len(pc_proteins)))
        self.set_curr_nodes_list()

    def Load_network_from_GS(self, genes):
        print('Gene Set, length:', len(genes))
        pc_proteins = [v.index for v in self.init_pa.graph.vs() if v['label'] in genes]
        pc_subgraph = self.init_pa.graph.induced_subgraph(pc_proteins)
        self.current_game_network = pc_subgraph
        self.current_game_network_filt = pc_subgraph
        self.CURR_SUBNET_SIZE.set(str(len(pc_proteins)))
        self.set_curr_nodes_list()

    def set_Ref_filter_value_sele(self, val):
        self.FILTER_REF_NUM_sele = float(val)
        try:
            sele = np.where(np.array(self.current_game_network.vs()[self.FILTER_TYPE.get()]) > self.FILTER_REF_NUM_sele)[0]
            cur_net_filt = self.current_game_network.induced_subgraph(sele)
        except Exception as e:
            print(e)
            print('Something Went Wrong')
        self.current_game_network_filt = cur_net_filt
        self.CURR_SUBNET_SIZE.set(str(len(self.current_game_network_filt.vs())))
        self.set_curr_nodes_list()

    def set_Ref_filter_value_init(self, val):
        self.FILTER_REF_NUM_init = float(val)
        try:
            sele_all = np.where(np.array(self.init_pa.graph.vs()[self.FILTER_TYPE.get()]) > self.FILTER_REF_NUM_init)[0]
            cur_init_net_filt = self.init_pa.graph.induced_subgraph(sele_all)
        except Exception as e:
            print(e)
            print('Something Went Wrong')
        self.init_network_filt = cur_init_net_filt
        self.CURR_ININET_SIZE.set(str(len(self.init_network_filt.vs())))


    def make_best_and_worst(self):
        df = self.STATS.stats_proteins.copy()
        df = df.drop(df[df.sum(axis=1) == 0].index)

        fig, ax = plt.subplots(1, 2, figsize=(9, 2))
        s1 = self.make_barplot(df, ax[0], False)
        s2 = self.make_barplot(df, ax[1], True)
        fig_x, fig_y = -470, -455
        fig_photo = draw_figure(self.canvas, fig, loc=(fig_x, fig_y))
        self.canvas.draw_idle()
        fig_w, fig_h = fig_photo.width(), fig_photo.height()

    def make_barplot(self, df, ax, ascending=False):
        df2 = df.sort_values(by='Score', ascending=ascending).iloc[:15]
        g = sns.barplot(x=df2.index, y=df2['Score'], data=df2, ax=ax)
        g.set_xticklabels(labels=df2.index, rotation=-45)
        g.tick_params(axis='both', which='major', labelsize=9.5)
        plt.gcf().subplots_adjust(bottom=0.33)
        lh = 'Highest'
        if ascending: lh = 'Lowest'
        g.set_title(lh + ' Score Known Proteins')
        return g

    def Show_uniprot2(self, ):
        """Depricated!"""
        print('HEY!')


