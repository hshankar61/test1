import empro
import sptk as sp
import matplotlib.pyplot as plt
import matplotlib
#matplotlib.use("Agg")

snpfile_name = r".\src\Connector_Model_KP4.s16p"
S1 = sp.sptk(snpfile_name)  # read in the S-parameter file
S1.set_port_config_se(config=0)  # specify the single ended port configuration
fig , ax = S1.plot_insertion_loss([1, 5])  # insertion loss for ports 1 and 2
plt.savefig(r".\src\InsertionLoss_1and5.pdf")  # save plot to file
# create a list of all ports on left side
port_list = list(range(1, S1.get_num_ports() + 1, 2))
fig , ax = S1.plot_insertion_loss(port_list)  # plot and title
plt.savefig(r".\src\InsertionLoss_allPorts.pdf")  # save plot to file
plt.show()
