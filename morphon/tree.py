import copy


class Node:
    def __init__(self, ident, value=None):
        self.ident = ident
        self.value = value
        self.parent = None
        self.children = []

    def __str__(self):
        return '<%s: %s %s %s>' % (str(self.ident), str(self.value),
            str(self.parent), str(self.children))

    def add_child(self, ident):
        if ident not in self.children:
            self.children.append(ident)

    def del_child(self, ident):
        if ident in self.children:
            self.children.remove(ident)


class Tree:
    def __init__(self):
        self.nodes = {}

    def __str__(self):
        s = '<Tree:\n'
        for ident in self.traverse():
            s += ' ' + str(self.nodes[ident]) + '\n'
        s += '>'
        return s

    def is_root(self, ident):
        return self.nodes[ident].parent is None

    def is_leaf(self, ident):
        return len(self.nodes[ident].children) == 0

    def is_fork(self, ident):
        return len(self.nodes[ident].children) > 1

    def is_bifurcation(self, ident):
        return len(self.nodes[ident].children) == 2

    def value(self, ident):
        return self.nodes[ident].value

    def parent(self, ident):
        return self.nodes[ident].parent

    def children(self, ident):
        return self.nodes[ident].children

    def degree(self, ident=None):
        leaves = filter(self.is_leaf, [item for item in self.traverse(ident)])
        return sum(1 for item in leaves)

    def order(self, ident):
        forks = filter(self.is_fork, [item for item in self.traverse(ident, reverse=True)])
        return sum(1 for item in forks)

    def root(self, ident=None):
        if len(self.nodes) == 0:
            return None
        ident = ident if ident is not None else self.nodes.iterkeys().next()
        path = [item for item in self.traverse(ident, reverse=True)]
        return path[-1]

    def traverse(self, ident=None, depth=True, reverse=False):
        if ident is None:
            ident = self.root()
        if ident in self.nodes:
            yield ident
            if reverse:
                parent = self.nodes[ident].parent
                while parent is not None:
                    yield parent
                    parent = self.nodes[parent].parent
            else:
                queue = self.nodes[ident].children
                while queue:
                    item = queue[0]
                    yield item
                    expansion = self.nodes[item].children
                    if depth:
                        queue = expansion + queue[1:]
                    else:
                        queue = queue[1:] + expansion

    def branch(self, ident=None, reverse=False):
        b = []
        item = ident if ident is not None else self.root()
        if item in self.nodes:
            b.append(item)
            if reverse:
                item = self.nodes[item].parent
                while item is not None and not self.is_fork(item):
                    b.insert(0, item)
                    item = self.nodes[item].parent
            else:
                while not self.is_leaf(item) and not self.is_fork(item):
                    item = self.nodes[item].children[0]
                    b.append(item)
        return b

    def branches(self, ident=None, depth=True, reverse=False):
        item = ident if ident is not None else self.root()
        if item in self.nodes:
            sec = self.branch(item, reverse=reverse)
            yield sec
            end = sec[0] if reverse else sec[-1]
            children = self.nodes[end].children
            parent = [self.nodes[end].parent] if not self.is_root(end) else []
            queue = parent if reverse else children
            while queue:
                item = queue[0]
                sec = self.branch(item, reverse=reverse)
                yield sec
                end = sec[0] if reverse else sec[-1]
                children = self.nodes[end].children
                parent = [self.nodes[end].parent] if not self.is_root(end) else []
                expansion = parent if reverse else children
                if depth:
                    queue = expansion + queue[1:]
                else:
                    queue = queue[1:] + expansion

    def add(self, ident=None, value=None, parent=None):
        if ident is None or ident < 1:
            if len(self.nodes) == 0:
                ident = 1
            else:
                ident = max(self.nodes.keys()) + 1
        if ident not in self.nodes:
            self.nodes[ident] = Node(ident, value)
            if parent is not None and parent in self.nodes:
                self.nodes[ident].parent = parent
                self.nodes[parent].add_child(ident)
        return ident

    def insert(self, ident, value=None):
        parent = self.nodes[ident].parent
        inset = max(self.nodes.keys()) + 1
        self.nodes[inset] = Node(inset, value)
        self.nodes[inset].children = [ident]
        self.nodes[inset].parent = parent
        if parent is not None:
            self.nodes[parent].del_child(ident)
            self.nodes[parent].add_child(inset)

    def remove(self, ident):
        parent = self.nodes[ident].parent
        children = self.nodes[ident].children
        if parent is not None:
            self.nodes[parent].del_child(ident)
            self.nodes[parent].children.extend(children)
        if children:
            for child in children:
                self.nodes[child].parent = parent

    def copy(self, ident=None):
        if ident is None:
            ident = self.root()
        tree = Tree()
        for item in self.traverse(ident):
            tree.nodes[item] = copy.deepcopy(self.nodes[item])
        tree.nodes[ident].parent = None
        return tree

    def prune(self, ident):
        if ident in self.nodes:
            parent = self.nodes[ident].parent
            if parent is not None:
                self.nodes[parent].del_child(ident)
            items = [item for item in self.traverse(ident)]
            for item in items:
                del self.nodes[item]

    def renumber(self, shift=0, continuous=False):
	if shift:
            for item in self.nodes:
                node = self.nodes[item]
                node.ident += shift
                if node.parent is not None:
                    node.parent += shift
                for child in node.children[:]:
                    node.children.append(child + shift)
                    node.children.remove(child)
            for item in self.nodes.keys():
                self.nodes[item+shift] = self.nodes[item]
                del self.nodes[item]
        elif continuous:
            idents = dict()
	    for ident, item in enumerate(self.traverse(), start=1):
	        idents[item] = ident
            tree = Tree()
	    for ident, item in enumerate(self.traverse(), start=1):
                tree.nodes[ident] = copy.deepcopy(self.nodes[item])
                tree.nodes[ident].parent = idents[self.nodes[item].parent] if self.nodes[item].parent is not None else None
		for c, child in enumerate(self.nodes[item].children):
                    tree.nodes[ident].children[c] = idents[child]
            self.nodes = tree.nodes

    def graft(self, ident, tree):
        stem = copy.deepcopy(tree)
        shift = max(self.nodes.keys())
        stem.renumber(shift)
        for item in stem.nodes:
            self.nodes[item] = stem.nodes[item]
        root = stem.root()
        self.nodes[ident].add_child(root)
        self.nodes[root].parent = ident
