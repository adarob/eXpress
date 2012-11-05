import heapq

# need to divide by average length
# output format: '# of leaves, # of nodes' followed by edges for each tree on line separated by semi-colon

class SimilarityMatrix:
  def __init__(self, n):
    self.M = {}
    self.partners = []
    self.heap = []
    self.deleted = set([])
    self.n = n
    for i in xrange(n):
      self.partners.append(set([]))

  #assumes i < j
  def increment(self, i, j, amt=1.):
    assert (i != j)
    pair = (i,j)
    if pair in self.M:
      self.M[pair] += amt
    else:
      self.partners[i].add(j)
      self.partners[j].add(i)
      self.M[pair] = amt
  
  def init_heap(self):
    self.heap = []
    for pair, sim in self.M.iteritems():
      heapq.heappush(self.heap, (-sim, pair))

  #assumes i < j
  def merge(self, i, j):
    del self.M[(i,j)]
    self.deleted.add(i)
    self.deleted.add(j)
    self.partners.append((self.partners[i] | self.partners[j]) - self.deleted)
    for k in (self.partners[i] & self.partners[j]) - self.deleted:
      i_key = (min(i,k), max(i,k))
      j_key = (min(j,k), max(j,k))
      n_key = (k, self.n)
      self.M[n_key] = (self.M[i_key] + self.M[j_key]) / 2.
      del self.M[i_key]
      del self.M[j_key]
      heapq.heappush(self.heap, (-self.M[n_key], n_key))
    for k in (self.partners[i] - self.partners[j]) - self.deleted:
      i_key = (min(i,k), max(i,k))
      n_key = (k, self.n)
      self.M[n_key] = self.M[i_key]
      del self.M[i_key]
      heapq.heappush(self.heap, (-self.M[n_key], n_key))
    for k in (self.partners[j] - self.partners[i]) - self.deleted:
      j_key = (min(j,k), max(j,k))
      n_key = (k, self.n)
      self.M[n_key] = self.M[j_key]
      del self.M[j_key]
      heapq.heappush(self.heap, (-self.M[n_key], n_key))
    self.n += 1

  def get_max(self):
    while True:
      if len(self.heap) == 0:
        return (None, None)
      val, pair = heapq.heappop(self.heap)
      if not (pair[0] in self.deleted or pair[1] in self.deleted):
        return pair
        
#initialize matrix
def initialize_matrix(infile):
  import pysam
  samfile = pysam.Samfile(infile)
  M = SimilarityMatrix(samfile.nreferences)
  lengths = samfile.lengths
  curr_qname = ''
  curr_targets = set([]) 
  n = 0
  for read in samfile:
    if read.qname != curr_qname:
      n += 1
      if n % 100000 == 0:
        print n
      if len(curr_targets) > 1:
        curr_targets = list(curr_targets)
        curr_targets.sort()
        for i in xrange(len(curr_targets)):
          x = curr_targets[i]
          for j in xrange(i+1, len(curr_targets)):
            y = curr_targets[j]
            M.increment(x,y, 2./(lengths[x] + lengths[y]))
      curr_targets = set([])
      curr_qname = read.qname

    if read.is_unmapped:
      continue

    if not read.is_paired:
      curr_targets.add(read.tid)
      continue

    if read.is_proper_pair and read.is_reverse:
      curr_targets.add(read.tid)
      continue
  return M
  
class Tree:
  def __init__(self, id, *children):
    self.id = id
    self.parent = None
    self.set_children(*children)

  def is_leaf(self):
    return len(self.children) == 0
  
  def set_children(self, *children):
    for child in children:
      child.parent = self
    self.children = children

  def leaves(self):
    if self.is_leaf():
      return [self.id]
    
    leaves = []
    for child in self.children:
      leaves += child.leaves()
    return leaves

  def to_string(self):
    s = ''
    for child in self.children:
      s += child.to_string()
      s += '%d,%d; ' % (self.id, child.id)
 
    return s
    
class Forest:
  def __init__(self, num_leaves):
    self.nodes = []
    for i in xrange(num_leaves):
      self.nodes.append(Tree(i))
  
  def add_node(self, id, *child_ids):
    self.nodes.append(Tree(id, *(self.nodes[child_id] for child_id in child_ids)))
  
  def insert_node(self, node):
    self.nodes.append(node)
  
  def to_string(self):
    s = ''
    num_leaves = 0
    for node in self.nodes:
      if node.parent == None and len(node.children) > 0:
        s += node.to_string() + '\n'
      num_leaves += node.is_leaf()
    s = '%d, %d\n' % (num_leaves, len(self.nodes)) + s
    return s
  
  def __len__(self):
    return len(self.nodes)

# destructive for M
def build_forest(M):
  F = Forest(M.n)
  M.init_heap()
  while True:
    i,j = M.get_max()
    if i == None:
      return F
    F.add_node(M.n, i, j)
    M.merge(i,j)
  return F
  
def build_forest2(infile):
  import pysam
  samfile = pysam.Samfile(infile)
  
  trees = [i for i in xrange(samfile.nreferences)]
  F = Forest(samfile.nreferences)
  
  curr_qname = ''
  trees_to_merge = set([]) 
  leaf_ids = set([])
  n = 0
  for read in samfile:
    if read.qname != curr_qname:
      n += 1
      if n % 100000 == 0:
        print n
      if len(trees_to_merge) > 1:
        id = len(F)
        F.add_node(id, *trees_to_merge)
        for leaf in F.nodes[-1].leaves():
          trees[leaf] = id 
      elif len(trees_to_merge) == 1:
        parents = set(F.nodes[leaf_id].parent for leaf_id in leaf_ids) 
        if len(parents) = 1:
          parent = parents.pop()
          if len(leaf_ids) < len(parent.children):
            u = Tree(len(F), *(F.nodes[leaf_id] for leaf_id in leaf_ids))
            F.insert_node(u)
            new_children = [child for child in p.children if (not child.id in leaf_ids)] + [u]
            parent.set_children(*new_children)
      trees_to_merge = set([])
      leaf_ids = set([])
      curr_qname = read.qname

    if read.is_unmapped:
      continue

    if (not read.is_paired) or (read.is_proper_pair and read.is_reverse):
      trees_to_merge.add(trees[read.tid])
      leaf_ids.add(read.tid)
      continue

  return F

#M = initialize_matrix('/home/adarob/experiments/express/simulation/hg19_ucsc_err_flip/hits.bam')
#F = build_forest(M)
F = build_forest2('/home/adarob/experiments/express/simulation/hg19_ucsc_err_flip/hits.bam')
file('forest_flip2.out','w').write(F.to_string())
