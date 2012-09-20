import heapq

# need to divide by average length
# output format: '# of leaves, # of nodes' followed by edges for each tree on line separated by semi-colon

class SimilarityMatrix:
  def __init__(self, n):
    self.M = {}
    self.partners = {}
    self.heap = []
    self.n = n
  
  #assumes i < j
  def increment(self, i, j, amt=1):
    pair = (i,j)
    self.M[pair] = self.M.get(pair, 0) + amt
  
  def init_heap(self):
    self.heap = []
    for pair, sim in self.M.iteritems():
      heapq.heappush(self.heap, (-sim, pair))

  #assumes i < j
  def merge(self, i, j):
    del self.M[(i,j)]
    for k in xrange(self.n):
      if i != k and j != k:
        i_key = (min(i,k), max(i,k))
        j_key = (min(j,k), max(j,k))
        n_key = (k, self.n)
        self.M[n_key] = (self.M[i_key] + self.M[j_key]) / 2.
        del self.M[i_key]
        del self.M[j_key]
        heapq.heappush(self.heap, (-self.M[n_key], n_key))
        self.n += 1
        
  def find_max(self):
    while True:
      if len(self.heap) == 0:
        return (None, None)
      val, pair = heapq.heappop(self.heap)
      if pair in self.M:
        return pair
        
#initialize matrix
def initialize_matrix(infile):
  import pysam
  samfile = pysam.Samfile(infile, 'r')
  M = SimilarityMatrix(samfile.nreferences)
  curr_qname = ''
  curr_targets = [] 
  for read in samfile:
    if read.qname != curr_qname:
      if len(curr_targets) > 1:
        curr_targets.sort()
        k = 0
        for i in curr_targets:
          k += 1
          for j in curr_targets[k:]:
            M.increment(i,j, 2./(samfile.lengths[i] + samfile.lengths[j]))
        curr_targets = []
    curr_qname = read.qname

    if read.is_unmapped:
      continue

    if not read.is_paired:
      curr_targets.append(read.tid)
      continue

    if read.is_proper_pair and read.is_reverse:
      curr_targets.append(read.tid)
      continue
  return M
  
class Tree:
  def __init__(self, id, *children):
    self.id = id
    self.children = children
    self.parent = None
    for child in children:
      child.parent = self

  def is_leaf(self):
    return len(self.children) == 0

  def to_string(self):
    s = ''
    for child in self.children:
      s += child.to_string()
      s += '%d,%d; ' % (self.id, child.id)
 
    return s
    
class Forest:
  def __init__(self):
    self.nodes = {}
  
  def add_node(self, id, *child_ids):
    self.nodes[id] = Tree(id, *(self.nodes[child_id] for child_id in child_ids))
  
  def to_string(self):
    s = ''
    for node in self.nodes.values():
      if node.parent == None:
        s += node.to_string() + '\n'
      num_leaves += node.is_leaf()
    s = '%d, %d' % (num_leaves, len(self.nodes)) + s
    return s

# destructive for C
def build_forest(C):
  F = Forest()
  while True:
    i,j = C.find_min()
    if i == None:
      return F
    F.add_node(C.n, i, j)
    C.merge(i,j)
    