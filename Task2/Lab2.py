from ete3 import Tree

t=Tree("spike_tree.nwk")
camel=t.search_nodes(name="MN514967.1_4902-8458.")[0]
t.set_outgroup(camel)
t.write(format=1,outfile="spike_tree_rooted.nwk")
print(t)