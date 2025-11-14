from ete3 import Tree

t=Tree("spike_tree_hosts.nwk")
camel=t.search_nodes(name="MN514967.1_Camelus_dromedarius_Camel")[0]
t.set_outgroup(camel)
t.write(format=1,outfile="spike_tree_hosts_rooted.nwk")
print(t)