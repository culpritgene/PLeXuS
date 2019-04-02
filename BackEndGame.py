from imports import *
from utils import *

class labirinth_game():
    def __init__(self, parent, path_length, mode, start_node):
        self.parent = parent

        self.path_length = path_length  ### number of steps to take
        self.mode = mode  ### What modes can we use? Directed / Different Granularities / Going in circles / ???
        self.start_node = start_node
        self.path_variants_ids = []
        self.path_variants = []  ### save nodes presented to user on current step
        self.path_corrects = []  ### save correctly chosen nodes in the list
        self.path_chosen = []  ### user decision

        self.game_on = True
        self.lives = self.parent.TOTAL_LIVES.get()  #### move to config subclass
        self.current_node = start_node

        self.PROT_INFO = {}
        self.MASKED_LABELS = []

    def get_index_in_subnetwork(self, gs):
        return self.parent.current_game_network_filt.vs()['label'].index(gs)

    def get_index_in_initnetwork(self, gs):
        return self.parent.init_network_filt.vs()['label'].index(gs)

    def start_game_12(self):
        ret2 = self.start_game_01()
        if ret2 != 0:
            self.start_game_02()

    def start_game_01(self):
        ret1 = self.generate_True_Path()
        print(ret1)
        if ret1 != 0:
            self.generate_wrong_vars()
            self.make_variants_gs()
            if self.parent.LABIRINTH_CONFIG.CONNECTGRAPH.get():
                msg = 'Making True Subgraph'
                print(msg)
                self.make_game_subgraph()
            elif not self.parent.LABIRINTH_CONFIG.CONNECTGRAPH.get():
                self.make_game_newgraph()
                self.make_fake_edges()
            self.make_game_subgraph_layout(_type=self.parent.LAYOUT_type.get())
            #             print('Loading hints from Uniprot...')
            #             self.load_uniprot_data()
            #             print(f'Loaded for {len(self.game_subgraph.vs())} proteins...')
            self.make_MASKED_LABELS()
            self.PROT_INFO.update({k: self.filter_ALL(k) for k in self.game_subgraph.vs()['label']})
            self.optimize_game_graph_layout()
            return 1
        else:
            print('Insufficient graph piece left')
            return 0

    def estimate_paths(self, cuid, max_l):
        trps = np.array(
            self.parent.current_game_network_filt.get_shortest_paths(cuid))  ### get farthest points in current subgraph
        trps_lens = np.array([len(x) for x in trps])
        res = np.where(trps_lens >= max_l)[0]
        if len(res) > 0:
            return (trps[np.random.choice(res)][:max_l + 1], 0)
        else:
            pp = np.random.choice(trps, p=(trps_lens ** 4) / sum(trps_lens ** 4))
            return (pp, max_l - len(pp))

    def generate_True_Path(self):
        max_l = self.parent.LABIRINTH_CONFIG.MAX_LENGTH.get()
        cuid = self.get_index_in_subnetwork(self.current_node)
        for i in range(12):
            trp, max_l2 = self.estimate_paths(cuid, max_l)
            if max_l2 == 0:
                self.true_path_ids = trp
                self.true_path = self.parent.current_game_network_filt.vs()[trp]['label']
                return 1
            else:
                cuid2 = trp[-1]
                for ii in range(10):
                    trp2, max_l3 = self.estimate_paths(cuid2, max_l2)
                    if max_l3 <= 0 and len(list(set(trp).intersection(set(trp2)))) == 1:
                        trp3 = list(trp) + list(trp2)[1:]
                        self.true_path_ids = trp3[:max_l]
                        self.true_path = self.parent.current_game_network_filt.vs()[trp3]['label']
                        return 2
        return 0

    def generate_wrong_vars(self, ):
        assert self.true_path, print('Firstly we should generate a track of true gene-symbols, than add false names')
        self.path_variants_ids.append((self.get_index_in_initnetwork(self.true_path[0])))
        for gs in self.true_path[1:]:
            ind = self.get_index_in_initnetwork(gs)
            self.all_paths = self.parent.init_network_filt.get_shortest_paths(
                ind)  ## spectrum of paths from true node to any other node
            all_distances = np.array([len(x) for x in self.all_paths])
            possible_randoms = np.where(all_distances == 0)[0]

            ## first we try to include only disconnected nodes (for now)
            if len(possible_randoms) > self.parent.LABIRINTH_CONFIG.BRANCHING_NUM.get() - 1:
                print(ind, 'all unconnected partners')
                pass
            else:
                possible_randoms = list(possible_randoms)
                possible_randoms += list(np.where(all_distances > self.parent.CONTRASTIVE_STEP)[0])
            sele = list(
                np.random.choice(possible_randoms, self.parent.LABIRINTH_CONFIG.BRANCHING_NUM.get() - 1, replace=False))
            self.path_variants_ids.append(tuple(sele + [ind]))

    def make_variants_gs(self, ):
        variants = []
        for sublist in self.path_variants_ids:
            if not isinstance(sublist, int):
                new_sublist = [self.parent.init_network_filt.vs()[node_id]['label'] for node_id in sublist]
            else:
                new_sublist = tuple([self.parent.init_network_filt.vs()[sublist]['label']])
            variants.append(tuple(new_sublist))
        self.path_variants = variants

    def unpack_variant_ids(self):
        unpacked = [node for sublist in self.path_variants[1:] for node in sublist]
        unpacked.insert(0, self.path_variants[0][0])
        unpacked_ids = [self.parent.init_network_filt.vs()['label'].index(x) for x in unpacked]
        return unpacked_ids

    def make_game_subgraph(self, ):
        unpacked_ids = self.unpack_variant_ids()
        self.game_subgraph = self.parent.init_network_filt.induced_subgraph(unpacked_ids, implementation='create_from_scratch')

    def make_game_newgraph(self, ):
        unpacked_ids = self.unpack_variant_ids()
        self.game_subgraph = igraph.Graph()
        for vs in self.parent.init_network_filt.vs()[unpacked_ids]:
            self.game_subgraph.add_vertex(**vs.attributes())

    def add_layers_info(self, ):
        layers = list(np.repeat(np.arange(1, len(self.path_variants[1:]) + 1), self.parent.LABIRINTH_CONFIG.BRANCHING_NUM.get()))
        mask = np.insert(layers, 0, 0)
        return layers

    def make_MASKED_LABELS(self, ):
        ### depending on self.parent.MASKING_AREA ...
        if self.parent.MASK_VOL.get()=='subgraph':
           self.MASKED_LABELS = self.game_subgraph.vs()['label']
        elif self.parent.MASK_VOL.get() == 'neighbors':
            ids = [v.index for v in self.parent.init_network_filt.vs() if v['name'] in self.game_subgraph.vs()['name']]
            ids2 = list(np.unique(np.concatenate(self.parent.init_network_filt.neighborhood(ids, 1))))
            self.MASKED_LABELS = [v['label'] for v in self.parent.init_network_filt.vs() if v.index in ids2]
        elif self.parent.MASK_VOL.get() == 'all':
            ids = [v.index for v in self.parent.init_pa.graph.vs() if v['name'] in self.game_subgraph.vs()['name']]
            ids2 = list(np.unique(np.concatenate(self.parent.init_pa.graph.neighborhood(ids, 2))))
            self.MASKED_LABELS = [v['label'] for v in self.parent.init_pa.graph.vs() if v.index in ids2]
        else:
            pass

    def load_uniprot_data(self, ):
        filter_general = lambda x: re.sub(r'\([^()]*\)', '', str(x))
        name_list = self.game_subgraph.vs()['name']
        all_data = []
        for name in name_list:
            print(name)
            data = get_uniprot_data(name).annotations['comment_function']
            data = filter_general(data)
            all_data.append(data)
        self.game_subgraph.vs()['Uniprot_func'] = all_data

    def make_fake_edges(self):
        for s0, s1 in zip(self.true_path, self.path_variants[1:]):
            core_id = self.game_subgraph.vs()['label'].index(s0)
            for s11 in s1:
                next_id = self.game_subgraph.vs()['label'].index(s11)
                if s11 in self.true_path:
                    corr = 'Correct'
                else:
                    corr = 'False'
                print('New-edge:', s0, '--->', s11, corr)
                self.game_subgraph.add_edge(core_id, next_id, **{'Fair': corr})

    def make_game_subgraph_layout(self, _type='reing'):
        if _type == 'reing':
            self.lay = self.game_subgraph.layout_fruchterman_reingold(repulserad=self.game_subgraph.vcount() ** 2,
                        maxiter=1200, area=self.game_subgraph.vcount() ** 2, coolexp=4, )
        elif _type == 'sugi':
            self.lay = self.game_subgraph.layout_sugiyama(self.add_layers_info(), hgap=2)
        elif _type == 'trees':
            self.lay = self.game_subgraph.layout_reingold_tilford(mode='OUT', )
        elif _type == 'circle':
            self.lay = self.game_subgraph.layout_reingold_tilford_circular(root=[6], mode='OUT',)
        elif _type == 'kawaii':
            self.lay = self.game_subgraph.layout_kamada_kawai(sigma=self.game_subgraph.vcount()/15, initemp=20,
                                                coolexp=0.96, kkconst=self.game_subgraph.vcount() ** 3)


    def optimize_game_graph_layout(self):
        bbox = igraph.BoundingBox(0, 0, self.parent.CANVAS_WIDTH / 2, self.parent.CANVAS_HEIGHT / 2)
        layout = self.lay
        layout.fit_into(bbox)
        layout.center(0, 0)
        edges = []
        for e in self.game_subgraph.es():
            edges.append(e.tuple)
        self.edges = np.array(edges)
        self.coords = np.array(layout.coords)

    def filter_ALL(self, lab):
        filter_general = lambda x: re.sub(r'\([^()]*\)', '', str(x))
        # sub_biological = lambda x,y: re.sub(r'{}'.format(x), '***', str(y))

        rep = dict((re.escape(k), '***') for k in self.MASKED_LABELS)
        try:
            del rep[lab]
        except:
            pass
        ### insert types of filtering here!
        pattern = re.compile("|".join(rep.keys()))

        ind = self.game_subgraph.vs()['label'].index(lab)
        annot = self.game_subgraph.vs()[ind]['Uniprot_func']
        annot = filter_general(annot)
        annot = pattern.sub(lambda m: rep[re.escape(m.group(0))], annot)
        return annot

    def force_decision(self, variants, true):
        decision = input(f'Please select adjacent to {self.current_node} protein among {variants}: ')
        self.path_chosen.append(variants[int(decision)])
        try:
            if int(decision) == variants.index(true):
                print('Correct! next step...')
                return True
            else:
                self.lives = self.lives - 1
                if self.lives > 0:
                    print(f'Sorry, this time wrong. Lets try again! ({self.lives} lives left)')
                    return False
                else:
                    print(f'Sorry, game over.')
                    self.game_on = False
                    return False
        except TypeError:
            print(f'Please, input integer from 0 to {len(variants)}')
