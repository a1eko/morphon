from builtins import list, dict, filter, all
import copy


class Node:
    def __init__(self, ident, value=None, parent=None):
        self.ident = ident
        self.value = value
        self.parent = parent
        self.children = list()

    def __str__(self):
        return '{}: {} {} {}'.format(str(self.ident), str(self.value),
            str(self.parent), str(self.children))

    def add_child(self, ident):
        if ident not in self.children:
            self.children.append(ident)

    def del_child(self, ident):
        if ident in self.children:
            self.children.remove(ident)


class Nodes(object):
    def __init__(self):
        self.nodes = dict()

    def __str__(self):
        return str([str(self.nodes[item]) for item in sorted(self.nodes)])

    def clear(self):
        return self.nodes.clear()

    def parent(self, ident):
        return self.nodes[ident].parent

    def children(self, ident):
        return list(self.nodes[ident].children)

    def value(self, ident):
        return self.nodes[ident].value

    def is_root(self, ident):
        return self.parent(ident) is None

    def is_leaf(self, ident):
        return len(self.children(ident)) == 0

    def is_fork(self, ident):
        return len(self.children(ident)) > 1

    def is_bifurcation(self, ident):
        return len(self.children(ident)) == 2

    def is_inner(self, ident):
        return not (self.is_root(ident) or self.is_fork(ident) or self.is_leaf(ident))

    def is_isolated(self, ident):
        return self.is_root(ident) and self.is_leaf(ident)

    def size(self, ident=None):
        if ident:
            length = len(list(self.traverse(ident)))
        else:
            length = len(self.nodes)
        return length

    def root(self, ident=None):
        ident = ident if ident is not None else list(self.nodes.keys())[0]
        path = list(self.traverse(ident, reverse=True))
        return path[-1]

    def roots(self):
        return filter(self.is_root, self.nodes)

    def leaves(self, ident=None):
        if ident:
            items = filter(self.is_leaf, self.traverse(ident))
        else:
            items = filter(self.is_leaf, self.nodes)
        return items

    def forks(self, ident=None, reverse=False):
        if ident:
            items = filter(self.is_fork, self.traverse(ident, reverse=reverse))
        else:
            items = filter(self.is_fork, self.nodes)
        return items

    def bifurcations(self, ident=None, reverse=False):
        if ident:
            items = filter(self.is_bifurcation, self.traverse(ident, reverse=reverse))
        else:
            items = filter(self.is_bifurcation, self.nodes)
        return items

    def isolates(self):
        return self.select(self.is_isolated)

    def degree(self, ident):
        return sum(1 for item in self.leaves(ident))

    def order(self, ident):
        return sum(1 for item in self.forks(ident, reverse=True))

    def add(self, node):
        self.nodes[node.ident] = node

    def connect(self, ident, parent):
        self.nodes[ident].parent = parent
        if parent is not None:
            self.nodes[parent].add_child(ident)

    def disconnect(self, ident):
        parent = self.parent(ident)
        self.nodes[ident].parent = None
        if parent is not None:
            self.nodes[parent].del_child(ident)

    def remove(self, ident):
        parent = self.parent(ident)
        children = self.children(ident)
        self.disconnect(ident)
        for child in children:
            self.disconnect(child)
            self.connect(child, parent)
        del self.nodes[ident]

    def insert(self, ident, node):
        node.children.clear()
        parent = self.parent(ident)
        self.add(node)
        if parent is not None:
            self.disconnect(ident)
            self.connect(node.ident, parent)
        self.connect(ident, node.ident)

    def copy(self, ident=None):
        if ident:
            n = Nodes()
            for item in self.traverse(ident):
                n.nodes[item] = copy.deepcopy(self.nodes[item])
            n.nodes[ident].parent = None
        else:
            n = copy.deepcopy(self)
        return n

    def prune(self, ident):
        parent = self.parent(ident)
        self.disconnect(ident)
        for item in list(self.traverse(ident)):
            del self.nodes[item]

    def traverse(self, ident=None, depth=True, reverse=False):
        ident = ident if ident is not None else self.root()
        yield ident
        if reverse:
            parent = self.parent(ident)
            while parent is not None:
                yield parent
                parent = self.parent(parent)
        else:
            queue = self.children(ident)
            while queue:
                item = queue[0]
                yield item
                expansion = self.children(item)
                if depth:
                    queue = expansion + queue[1:]
                else:
                    queue = queue[1:] + expansion

    def section(self, ident, reverse=False):
        for item in self.traverse(ident, reverse=reverse):
            yield item
            term = item if not reverse else self.parent(item)
            if term and not self.is_inner(term):
                break

    def sections(self, ident=None, depth=True, reverse=False, with_parent=False):
        ident = ident if ident is not None else self.root()
        sec = list(self.section(ident, reverse=reverse))
        #
        def add_parent():
            if reverse:
                parent = self.parent(sec[-1])
                if parent is not None:
                    sec.append(parent)
            else:
                parent = self.parent(sec[0])
                if parent is not None:
                    sec.insert(0, parent)
        #
        if with_parent:
            add_parent()
        yield sec
        if reverse:
            item = sec[-1]
            parent = self.parent(item) if not with_parent else item
            while not self.is_root(item) and not self.is_root(parent):
                sec = list(self.section(parent, reverse=True))
                if with_parent:
                    add_parent()
                yield sec
                item = sec[-1]
                parent = self.parent(item) if not with_parent else item
        else:
            item = sec[-1]
            queue = self.children(item)
            while queue:
                item = queue[0]
                sec = list(self.section(item))
                if with_parent:
                    add_parent()
                yield sec
                item = sec[-1]
                expansion = self.children(item)
                if depth:
                    queue = expansion + queue[1:]
                else:
                    queue = queue[1:] + expansion

    def renumber(self, shift=0):
        n = Nodes()
        if shift:
            for (item, node) in self.nodes.items():
                n.nodes[item+shift] = copy.deepcopy(node)
                n.nodes[item+shift].children = list()
                for child in node.children:
                    n.nodes[item+shift].add_child(child+shift)
                n.nodes[item+shift].ident = node.ident + shift
                if node.parent is not None:
                    n.nodes[item+shift].parent = node.parent + shift
        else:
            idents = dict()
            for ident, item in enumerate(self.traverse(), start=1):
                idents[item] = ident
            for ident, item in enumerate(self.traverse(), start=1):
                n.nodes[ident] = copy.deepcopy(self.nodes[item])
                n.nodes[ident].parent = idents[self.parent(item)] if not self.is_root(item) else None
                for c, child in enumerate(self.children(item)):
                    n.nodes[ident].children[c] = idents[child]
                n.nodes[ident].ident = ident
        self.nodes = n.nodes

    def graft(self, ident, nodes):
        shift = max(self.nodes.keys())
        stem = copy.deepcopy(nodes)
        stem.renumber(shift)
        for item in stem.nodes:
            self.nodes[item] = stem.nodes[item]
        root = stem.root()
        self.nodes[ident].add_child(root)
        self.nodes[root].parent = ident

    def check(self):
        rep = dict()
        roots = list(self.roots())
        nodes = set(self.nodes[item].ident for item in self.nodes)
        parents = set(self.parent(item) for item in filter(lambda i: not self.is_root(i), self.nodes))
        children = list()
        for item in self.nodes:
            children.extend(self.children(item))
        children = set(children)
        rep['nonempty'] = self.size() > 0
        rep['single_root'] = len(roots) == 1
        rep['defined_ids'] = all(nodes)
        rep['unique_ids'] = len(nodes) == len(self.nodes)
        rep['valid_ids'] = nodes == set(self.nodes)
        rep['valid_parents'] = parents.issubset(self.nodes)
        rep['valid_children'] = children.issubset(self.nodes)
        rep['no_root_reference'] = children.isdisjoint(roots)
        rep['no_self_reference'] = all(self.parent(item) not in self.children(item)
            for item in self.nodes)
        rep['valid_links'] = all(item in self.children(self.parent(item))
            for item in self.nodes if self.parent(item) is not None)
        rep['traversable_data'] = set(self.nodes) == set(self.traverse())
        return rep

    def is_tree(self):
        return all(self.check().values())
