import os
import empro
import matplotlib.pyplot as plt
import numpy as np


class sptk:
    """
    Implement features similar to that found in the S-Parameter Tool Kit (sptk) in ADS
    """
    def __init__(self, snpfile):
        """Initialize the sptk class by reading in the S-Parameter file.
        The class can be in one of 2 states - single ended ('se') or differential
        ('diff'). They are set by calling the methods compute_se_data() or
        compute_diff_data() respectively. In the 'se' state, the S parameters
        are all single ended. In the 'diff' state, the S Parameters are mixed
         mode, differntial and common. To go into the 'se' state, call the method
        compute_se_data(). To go into the 'diff' state, call the method
        compute_diff_data(). Also associated with each of these states
        is a configuration, config_se and config_diff. These are set by calling
        the methods set_port_config_se() or set_port_config_diff() respectively.

        Args:
            snpfile (string): This is the name of the S-Parameter file

        Raises:
            ValueError: The number of ports must be even (since every port is connected
            to another port called its through port)
        """
        self.spar_file_name = snpfile  # save input file name
        # read in file and parse
        self._results = empro.analysis.CircuitResults(snpfile,
                                                      os.path.split(snpfile)[1])
        self.num_ports = self._results.numberOfPorts()  # number of single ended ports
        if self.num_ports % 2 == 1:
            raise ValueError("Number of ports must be even")
        self.freq = np.array(self._results.frequencies())  # array of frequencies
        self.config_se = 0
        self.config_diff = 0  # configuration for diff ports
        self.num_diff_ports = self.num_ports // 2  # num_ports must be even
        # The class can be in one of 2 states - single ended or differential
        # Thie state is set by a call to the functions compute_se_data() or
        # compute_diff_data(). By default state == "se"
        self.state = "se"

    def set_port_config_se(self, config=0):
        """ This function is used to specify how the S-Parameter ports are connected.
        If port i is connected to port j, then j is the through port of port i and
        i is the through port of port j. There also is a side, left or right,
        associated with each port.
        There are different ways of assigning through ports as in Sparameter Toolkit.
        Here is an example if there are 8 ports
        config = 0 implies the connection 1->2, 3->4, 5->6, 7->8,
        ports [1,3,5,7] are on the left side, [2,4,6,8] are on the right side
        config = 1 impiles the connection 1->5, 2->6, 3->7, 4->8,
        ports [1,2,3,4] are on the left side, [5,6,7,8] are on the right side
        config = 2 implies the connection 1->8, 2->7, 3->6, 4->5,
        ports [1,2,3,4] are on the left side, [8,7,6,5] are on the right side
        See the documentation for figures.
        Args:
            config (int, optional): Defaults to 0.

        Raises:
            ValueError: The parameter config must take one of the values [0, 1, 2],
            else error
        """
        if (config < 0) or (config > 2):
            raise ValueError("config parameter must be 0 or 1 or 2")
        self.config_se = config

    def set_port_config_diff(self, config=0):
        """This method combines single ended ports into differential ports.
        This combination of pairs of single ended ports is
        based on the single ended port configuration.
        The labeling of the differentail ports is controlled by the config parameter.
        See the documentation for figures.

        Args:
            config (int, optional):  Defaults to 0.

        Raises:
            ValueError:  The parameter config must take one of the values
            [0, 1, 2], else error
        """
        # Assumes that assign_through_ports() has been called already
        # and self.config is set
        if (config < 0) or (config > 2):
            raise ValueError("config parameter must be 0 or 1 or 2")
        self.config_diff = config

    def get_num_ports(self, state=None):
        """Returns the number of ports. Depending on the state of the class
        or the input parameter state, the number of single ended ports or
        the number of differential ports is returned

        Args:
            state (_type_, optional): Can take the values None, 'se' or 'diff'.
            If 'None', the self.state value is used.

        Raises:
            ValueError: The parameter state takes the values None, 'se' or 'diff'

        Returns:
            int: number of ports, single or differential
        """
        if state is None:
            state = self.state
        if state == "se":
            return self.num_ports
        elif state == "diff":
            return self.num_diff_ports
        else:
            raise ValueError("not a valid value for self.state")

    def compute_diff_data(self):
        """Computes the mixed mode S-parameters, diff and common mode and
        sets the state of the class as "diff"

        Raises:
            ValueError: If the state has been set to a value other than "se" or "diff",
            then error
        """
        if self.num_ports % 4 != 0:
            raise ValueError("Number of ports must be a multiple of \
                              4 for differential ports")
        diff_pairs = []
        # diff_pairs is a list of differential ports
        diff_pairs = [self.diff_to_se(k) for k in range(1, self.num_diff_ports + 1)]
        # Subtract 1 from dif_pairs because port numbers in _results start at 0
        diff_pairs = [[y - 1 for y in x] for x in diff_pairs]
        self._results.setupDifferentialPairs(diff_pairs)
        self.state = "diff"

    def compute_se_data(self):
        """Computes the single ended S-parameters and
        sets the state of the class as "se"
        """
        # Resets any ports set as differential to single ended
        # empty list sets Spars to single ended
        self._results.setupDifferentialPairs([])
        self.state = "se"

    def get_through_port(self, port_number, state=None):
        """Find the corresponding through port of port_number, i.e.,
        which port is it connected to on the other side

        Args:
            port_number (int): _description_
            state (None, 'se', 'diff', optional): Get the corresponding through port
            which depends on the state and the config value. Defaults to None.

        Raises:
            ValueError: Chceks that the port_number and state are valid values

        Returns:
            int: retruns the port that the input parameter port_number is connected to
        """
        # Find the corresponding through port of port_number, i.e., which port
        # is it connected to on the other side
        # Return values are the corresponding through port and the side it is on
        # side 0 is left side and side 1 is right side
        if state is None:
            state = self.state
        if state == 'se':
            num_ports = self.get_num_ports(state)
            config = self.config_se
        elif state == "diff":
            num_ports = self.get_num_ports(state)
            config = self.config_diff
        else:
            raise ValueError("Not in a valid state")

        if port_number < 1 or port_number > num_ports:
            raise ValueError("Not a valid port number")

        if config == 0:
            if (port_number % 2 == 0):  # port_number is on right, through port on left
                return port_number - 1, 0
            else:
                return port_number + 1, 1
        elif config == 1:
            if (port_number <= num_ports // 2):
                return port_number + num_ports // 2, 1
            else:
                return port_number - num_ports // 2, 0
        elif config == 2:
            if (port_number <= num_ports // 2):
                return num_ports - port_number + 1, 1
            else:
                return num_ports - port_number + 1, 0
        else:
            raise ValueError("self.config is not a valid value")

    def diff_to_se(self, diff_port_number):
        """Returns the 2 single ports which have been combined as a differential port.
        This depends on the config values set by the methods set_port_config_se()
        and set_port_config_diff()

        Args:
            diff_port_number (int): Differential port number

        Raises:
            ValueError: Checks is the diff_port_number and configuration
            values are valid


        Returns:
            int: a list of integers which are the single ended ports
            associated with the differential port
        """
        if diff_port_number < 1 or diff_port_number > self.num_diff_ports:
            raise ValueError("not a valid differential port number")
        num_diff_ports = self.num_diff_ports
        if (self.config_se == 0 and self.config_diff == 0):
            foo = 2 * diff_port_number
            if diff_port_number % 2 == 1:  # left side diff port
                return [foo - 1, foo + 1]
            else:  # right side diff port
                return [foo - 2, foo]
        elif (self.config_se == 0 and self.config_diff == 1):
            if (diff_port_number <= num_diff_ports // 2):  # left side diff port
                foo = (diff_port_number - 1) * 4 + 1
                return [foo, foo + 2]
            else:
                foo = 4 * diff_port_number - self.num_ports
                return [foo - 2, foo]
        elif (self.config_se == 0 and self.config_diff == 2):
            if (diff_port_number <= num_diff_ports // 2):  # left side diff port
                foo = 4 * (diff_port_number - 1)
                return [foo + 1, foo + 3]
            else:
                foo = 2 * self.num_ports - 4 * diff_port_number + 4
                return [foo - 2, foo]
        elif (self.config_se == 1 and self.config_diff == 0):
            if diff_port_number % 2 == 1:  # left side diff port
                return [diff_port_number, diff_port_number + 1]
            else:
                foo = num_diff_ports + diff_port_number
                return [foo - 1, foo]
        elif (self.config_se == 1 and self.config_diff == 1):
            foo = 2 * diff_port_number
            return [foo - 1, foo]
        elif (self.config_se == 1 and self.config_diff == 2):
            if (diff_port_number <= num_diff_ports // 2):  # left side diff port
                foo = 2 * diff_port_number
                return [foo - 1, foo]
            else:
                foo = 3 * num_diff_ports - 2 * diff_port_number
                return [foo + 1, foo + 2]
        elif (self.config_se == 2 and self.config_diff == 0):
            if diff_port_number % 2 == 1:  # left side diff port
                return [diff_port_number, diff_port_number + 1]
            else:
                foo = self.num_ports - diff_port_number
                return [foo + 2, foo + 1]
        elif (self.config_se == 2 and self.config_diff == 1):
            if (diff_port_number <= num_diff_ports // 2):  # left side diff port
                return [2 * diff_port_number - 1, 2 * diff_port_number]
            else:
                foo = self.num_ports - 2 * diff_port_number + num_diff_ports
                return [foo + 2, foo + 1]
        elif (self.config_se == 2 and self.config_diff == 2):
            if (diff_port_number <= num_diff_ports // 2):  # left side diff port
                return [2 * diff_port_number - 1, 2 * diff_port_number]
            else:
                foo = (diff_port_number - num_diff_ports // 2 - 1) * 2 + num_diff_ports
                return [foo + 2 , foo + 1]
        else:
            raise ValueError("config values are not valid")

    def get_ne_ports(self, port_number, state=None):
        """Find all the ports on the same side as port_number , i.e., near end ports.
        Depending on the configuration set by the methods set_port_config_se() and
        set_port_config_diff(), each port
        is assumed to be on the left or right side (see documentation for figures)

        Args:
            port_number (int): The port number with respect to which the
            near end ports are found
            state (None, 'se' or 'diff', optional): The configuration
            showing the location and connection of ports. Defaults to None.

        Raises:
            ValueError: If port_number or state are not valid, then error

        Returns:
            list of ints: a list of port numbers on the same side as port_number
        """
        #
        if state is None:
            state = self.state
        if state == 'se':
            num_ports = self.get_num_ports(state)
            config = self.config_se
        else:  # 'diff'
            num_ports = self.get_num_ports(state)
            config = self.config_diff
        if port_number < 1 or port_number > num_ports:
            raise ValueError("not a valid port number")

        if config == 0:
            if (port_number % 2 == 0):  # even port_numbers are on right side
                ne_ports = [x for x in range(2, num_ports + 1, 2)]
                ne_ports.remove(port_number)
            else:
                ne_ports = [x for x in range(1, num_ports, 2)]
                ne_ports.remove(port_number)
        elif config == 1:
            if (port_number <= self.num_ports // 2):
                ne_ports = [x for x in range(1, num_ports // 2 + 1, 1)]
                ne_ports.remove(port_number)
            else:
                ne_ports = [x for x in range(num_ports // 2 + 1, num_ports + 1, 1)]
                ne_ports.remove(port_number)
        elif config == 2:
            if (port_number <= num_ports // 2):
                ne_ports = [x for x in range(1, num_ports // 2 + 1, 1)]
                ne_ports.remove(port_number)
            else:
                ne_ports = [x for x in range(num_ports, num_ports // 2, -1)]
                ne_ports.remove(port_number)
        else:
            raise ValueError("self.config is not a valid value")
        return ne_ports

    def get_fe_ports(self, port_number, state=None):
        """Find all the ports on the oppsite side as port_number , i.e., far end ports.
        Depending on the configuration set by the methods set_port_config_se() and
        set_port_config_diff(), each port
        is assumed to be on the left or right side (see documentation for figures)

        Args:
            port_number (int): The port number with respect to which the far
            end ports are found
            state (None, 'se' or 'diff', optional): The configuration showing
            the location and connection of ports. Defaults to None.

        Raises:
            ValueError: If port_number or state are not valid, then error

        Returns:
           list of ints: a list of port numbers on the opposite side as port_number
        """
        if state is None:
            state = self.state
        # Find all the ports on the other side as port_number , i.e., far end ports
        if state == 'se':
            if port_number < 1 or port_number > self.num_ports:
                raise ValueError("not a valid port number")
            num_ports = self.num_ports
            config = self.config_se
        else:  # 'diff'
            if port_number < 1 or port_number > self.num_diff_ports:
                raise ValueError("not a valid diff port number")
            num_ports = self.num_ports // 2
            config = self.config_diff

        if config == 0:
            if (port_number % 2 == 0):  # even port_numbers are on right side
                ne_ports = [x for x in range(1, num_ports, 2)]
                ne_ports.remove(self.get_through_port(port_number)[0])
            else:
                ne_ports = [x for x in range(2, num_ports + 1, 2)]
                ne_ports.remove(self.get_through_port(port_number)[0])
        elif config == 1:
            if (port_number <= self.num_ports // 2):
                ne_ports = [x for x in range(num_ports // 2 + 1, num_ports + 1, 1)]
                ne_ports.remove(self.get_through_port(port_number)[0])
            else:
                ne_ports = [x for x in range(1, num_ports // 2 + 1, 1)]
                ne_ports.remove(self.get_through_port(port_number)[0])
        elif config == 2:
            if (port_number <= num_ports // 2):
                ne_ports = [x for x in range(num_ports, num_ports // 2, -1)]
                ne_ports.remove(self.get_through_port(port_number)[0])
            else:
                ne_ports = [x for x in range(1, num_ports // 2 + 1, 1)]
                ne_ports.remove(self.get_through_port(port_number)[0])
        else:
            raise ValueError("self.config is not a valid value")
        return ne_ports

    def plot(self, port_pairs, plot_title=""):
        """Plots the S parameters associated with a list of pair of ports

        Args:
            port_pairs (list of lists): [[p1, p2], [p3,p4]]
            plot_title (str, optional): Title string used in plot. Defaults to "".

        Returns:
            fig, ax: Matplotlib objects
        """
        fig, ax = plt.subplots(1, 1)
        plt.xlabel('F (Hz)')
        plt.ylabel('dB')
        plt.grid()
        plt.title(plot_title)
        for k in port_pairs:
            S_mag = np.array(self._results.Src(k[1] - 1, k[0] - 1, "ComplexMagnitude"))
            plt.semilogx(np.array(self._results.frequencies()), 20 * np.log10(S_mag),
                         label='S' + str(k[1]) + ',' + str(k[0]))
        plt.legend()
        ax = plt.gca()
        return fig, ax

    def plot_insertion_loss_se(self, port_list, plot_title='Insertion Loss'):
        """Plots the single ended insertion loss for a list of pair of ports.
        If this function is called and the class is in the 'diff' state,
        the state is changed to 'se'

        Args:
            port_pairs (list of lists): [[p1, p2], [p3,p4]]
            plot_title (str, optional): Title string used in plot.
            Defaults to "Insertion Loss".

        Raises:
            ValueError: Throws an error if an invalid port number is in port_list

        Returns:
           fig, ax: Matplotlib objects
        """
        if self.state == "diff":
            self.compute_se_data()
        fig, ax = plt.subplots(1, 1)
        plt.xlabel('F (Hz)')
        plt.ylabel('dB')
        plt.grid()
        plt.title(plot_title)
        for k in port_list:
            if (k < 1) or (k > self.get_num_ports()):
                raise ValueError("not a valid port number")
            p_through, p_through_side = self.get_through_port(k)
            Sp_through_k_mag = np.array(self._results.Src(p_through - 1,
                                        k - 1, "ComplexMagnitude"))
            plt.semilogx(np.array(self._results.frequencies()),
                         20 * np.log10(Sp_through_k_mag),
                         label='S' + str(p_through) + ',' + str(k))
        plt.legend()
        ax = plt.gca()
        return fig, ax

    def plot_insertion_loss_diff(self, diff_port_list, plot_diff=True,
                                 plot_mixed=False, plot_title="Insertion Loss Diff"):
        """Plots the mixed mode, diff and common,
        insertion loss fo a list of pair of ports.
        If this function is called and the class is in the 'se' state,
        the state is changed to 'diff'

        Args:
            diff_port_list (list of lists): [[p1, p2], [p3,p4]]
            plot_diff (True or False): plot the diff and common mode terms.
            Defaults to True.
            plot_mixed (True or False): plot the mode conversion terms, diff to
            common and common to diff. Defaults to False
            plot_title (str, optional): Title string used in plot.
            Defaults to "Insertion Loss Diff".

        Raises:
            ValueError: Throws an error if an invalid port number is in port_list

        Returns:
            fig, ax: Matplotlib objects
        """
        if self.state == "se":
            self.compute_differential_data()
        fig, ax = plt.subplots(11)
        plt.xlabel('F (Hz)')
        plt.ylabel('dB')
        plt.grid()
        plt.title(plot_title)
        for k in diff_port_list:
            if (k < 1) or (k > self.get_num_ports()):
                raise ValueError("not a valid port number")
            k1 = self.diff_to_se(k)  # returns a pair of ports [k11, k12]
            k11 = k1[0]
            k11_through, _ = self.get_through_port(k11, state="se")
            k12 = k1[1]
            k12_through, _ = self.get_through_port(k12, state="se")
            # At this point we have the diff pair
            # k11 --- k11_through
            # k12 --- k12_through
            if plot_diff is True:
                Sd2d1 = np.array(self._results.Src(k11_through - 1,
                                 k11 - 1, "ComplexMagnitude"))
                k_diff_through, _ = self.get_through_port(k)
                plt.semilogx(np.array(self._results.frequencies()),
                             20 * np.log10(Sd2d1),
                             label='Sd' + str(k_diff_through) + 'd' + str(k))
                Sc2c1 = np.array(self._results.Src(k12_through - 1,
                                 k12 - 1, "ComplexMagnitude"))
                plt.semilogx(np.array(self._results.frequencies()),
                             20 * np.log10(Sc2c1),
                             label='Sc' + str(k_diff_through) + 'c' + str(k))
            if plot_mixed is True:
                Sc2d1 = np.array(self._results.Src(k12_through - 1,
                                 k11 - 1, "ComplexMagnitude"))
                Sd2c1 = np.array(self._results.Src(k11_through - 1,
                                 k12 - 1, "ComplexMagnitude"))
                plt.semilogx(np.array(self._results.frequencies()),
                             20 * np.log10(Sc2d1),
                             label='Sc' + str(k_diff_through) + 'd' + str(k))
                plt.semilogx(np.array(self._results.frequencies()),
                             20 * np.log10(Sd2c1),
                             label='Sd' + str(k_diff_through) + 'c' + str(k))
        plt.legend()
        ax = plt.gca()
        return fig, ax

    def plot_insertion_loss(self, port_list, plot_title='Insertion Loss',
                            plot_diff=True, plot_mixed=False):
        """Plots the insertion loss for a list of pair of ports.
        Depending on the state, 'se' or 'diff', either the single ended or
        diff mode S-parameters are plotted

        Args:
            port_list (list of lists): [[p1, p2], [p3,p4]], all ints
            plot_title (str, optional): Title used for plot.
            Defaults to 'Insertion Loss'.
            plot_diff (bool, optional): Used only if state == 'diff' and
            then plots diff to diff and common to common mode components.
            Defaults to True.
            plot_mixed (bool, optional): Used only if state == 'diff' and
            then plots diff to common and common to diff mode components.
            Defaults to False.

        Raises:
            ValueError: Throws an error if an invalid port number is
            specified in port_list

        Returns:
            _type_: fig, ax: Matplotlib objects
        """
        if self.state == "se":
            fig, ax = self.plot_insertion_loss_se(port_list, plot_title=plot_title)
        elif self.state == "diff":
            fig, ax = self.plot_insertion_loss_diff(port_list, plot_title=plot_title,
                                                    plot_diff=plot_diff,
                                                    plot_mixed=plot_mixed)
        else:
            raise ValueError("State self.state has an invalid value")
        return fig, ax

    def plot_return_loss_se(self, port_list, plot_title='Return Loss'):
        """Plots the single ended return loss for a list of pair of ports.
        If this function is called and the class is in the 'diff' state,
        the state is changed to 'se'

        Args:
            port_pairs (list of lists): [[p1, p2], [p3,p4]]
            plot_title (str, optional): Title string used in plot.
            Defaults to "Return Loss".

        Raises:
            ValueError: Throws an error if an invalid port number is in port_list

        Returns:
           fig, ax: Matplotlib objects
        """
        if self.state == "diff":
            self.compute_se_data()
        fig, ax = plt.subplots(1, 1)
        plt.xlabel('F (Hz)')
        plt.ylabel('dB')
        plt.grid()
        plt.title(plot_title)
        for k in port_list:
            if (k <= 0) or (k > self.get_num_ports()):
                raise ValueError("not a valid port number")
            Sp_return_k_mag = np.array(self._results.Src(k - 1,
                                       k - 1, "ComplexMagnitude"))
            plt.semilogx(np.array(self._results.frequencies()),
                         20 * np.log10(Sp_return_k_mag),
                         label='S' + str(k) + ',' + str(k))
        plt.legend()
        return fig, ax

    def plot_return_loss_diff(self, port_list, plot_diff=True, plot_mixed=False,
                              plot_title='Return Loss Diff'):
        """Plots the mixed mode, diff and common, return loss
        for a list of pair of ports.
        If this function is called and the class is in the 'se' state,
        the state is changed to 'diff'

        Args:
            diff_port_list (list of lists): [[p1, p2], [p3,p4]]
            plot_diff (True or False): plot the diff and common mode terms.
            Defaults to True
            plot_mixed (True or False): plot the mode conversion terms,
            diff to common and common to diff. Defaults to False
            plot_title (str, optional): Title string used in plot.
            Defaults to "Return Loss Diff".

        Raises:
            ValueError: Throws an error if an invalid port number is in port_list

        Returns:
            fig, ax: Matplotlib objects
        """
        if self.state == "se":
            self.compute_differential_data()
        fig, ax = plt.subplots(1, 1)
        plt.xlabel('F (Hz)')
        plt.ylabel('dB')
        plt.grid()
        plt.title(plot_title)
        for k in port_list:
            if (k < 1) or (k > self.get_num_ports()):
                raise ValueError("not a valid port number")
            se_ports = self.diff_to_se(k)
            if plot_diff is True:
                Sp_return_k_mag = np.array(self._results.Src(se_ports[0] - 1,
                                           se_ports[0] - 1, "ComplexMagnitude"))
                plt.semilogx(np.array(self._results.frequencies()),
                             20 * np.log10(Sp_return_k_mag),
                             label='Sd' + str(k) + 'd' + str(k))
                Sp_return_k_mag = np.array(self._results.Src(se_ports[1] - 1,
                                           se_ports[1] - 1, "ComplexMagnitude"))
                plt.semilogx(np.array(self._results.frequencies()),
                             20 * np.log10(Sp_return_k_mag),
                             label='Sc' + str(k) + 'c' + str(k))
            if plot_mixed is True:
                Sp_return_k_mag = np.array(self._results.Src(se_ports[0] - 1,
                                           se_ports[1] - 1, "ComplexMagnitude"))
                plt.semilogx(np.array(self._results.frequencies()),
                             20 * np.log10(Sp_return_k_mag),
                             label='Sd' + str(k) + 'c' + str(k))
                Sp_return_k_mag = np.array(self._results.Src(se_ports[1] - 1,
                                           se_ports[0] - 1, "ComplexMagnitude"))
                plt.semilogx(np.array(self._results.frequencies()),
                             20 * np.log10(Sp_return_k_mag),
                             label='Sc' + str(k) + 'd' + str(k))

        plt.legend()
        return fig, ax

    def plot_return_loss(self, port_list, plot_title='Return Loss', plot_diff=True,
                         plot_mixed=False):
        """Plots the return loss associated with a list of pair of ports. \n
        Depending on the state, 'se' or 'diff', either the single ended or diff mode
        S-parameters are plotted

        Args:
            port_list (list of lists): [[p1, p2], [p3,p4]], all ints
            plot_title (str, optional): Title used for plot. Defaults to 'Return Loss'.
            plot_diff (bool, optional): Used only if state == 'diff' and
            then plots diff to diff and common to common mode components.
            Defaults to True.
            plot_mixed (bool, optional): Used only if state == 'diff' and
            then plots diff to common and common to diff mode components.
            Defaults to False.

        Raises:
            ValueError: Throws an error if an invalid port number is
            specified in port_list

        Returns:
            _type_: fig, ax: Matplotlib objects
        """
        if self.state == "se":
            fig, ax = self.plot_return_loss_se(port_list, plot_title=plot_title)
        elif self.state == "diff":
            fig, ax = self.plot_return_loss_diff(port_list, plot_title=plot_title,
                                                 plot_diff=plot_diff,
                                                 plot_mixed=plot_mixed)
        else:
            raise ValueError("State self.state has an invalid value")
        return fig, ax

    def plot_next_victim_se(self, port_list, plot_title='Next Victim'):
        """Plots the single ended next victim coupling for a list of pair of ports.
        If this function is called and the class is in the 'diff' state,
        the state is changed to 'se'

        Args:
            port_pairs (list of lists): [[p1, p2], [p3,p4]]
            plot_title (str, optional): Title string used in plot.
            Defaults to "Next Victim".

        Raises:
            ValueError: Throws an error if an invalid port number is in port_list

        Returns:
           fig, ax: Matplotlib objects
        """
        if self.state == "diff":
            self.compute_se_data()
        fig, ax = plt.subplots(1, 1)
        plt.xlabel('F (Hz)')
        plt.ylabel('dB')
        plt.grid()
        plt.title(plot_title)
        for k0 in port_list:
            if (k0 < 1) or (k0 > self.get_num_ports()):
                raise ValueError("not a valid port number")
            p_ne = self.get_ne_ports(k0)
            for k1 in p_ne:
                S_k0_k1 = np.array(self._results.Src(k0 - 1,
                                   k1 - 1, "ComplexMagnitude"))
                plt.semilogx(self._results.frequencies(),
                             20 * np.log10(S_k0_k1),
                             label='S' + str(k0) + ',' + str(k1))
        plt.legend()
        return fig, ax

    def plot_next_victim_diff(self, port_list, plot_diff=True, plot_mixed=False,
                              plot_title='Next Victim Diff'):
        """Plots the mixed mode, diff and common, next victim coupling
        for a list of pair of ports.
        If this function is called and the class is in the 'se' state,
        the state is changed to 'diff'

        Args:
            diff_port_list (list of lists): [[p1, p2], [p3,p4]]
            plot_diff (True or False): plot the diff and common mode terms.
            Defaults to True
            plot_mixed (True or False): plot the mode conversion terms,
            diff to common and common to diff. Defaults to False
            plot_title (str, optional): Title string used in plot.
            Defaults to "Next Victim Diff".

        Raises:
            ValueError: Throws an error if an invalid port number is in port_list

        Returns:
            fig, ax: Matplotlib objects
        """
        if self.state == "se":
            self.compute_differential_data()
        fig, ax = plt.subplots(1, 1)
        plt.xlabel('F (Hz)')
        plt.ylabel('dB')
        plt.grid()
        plt.title(plot_title)
        for k0 in port_list:
            if (k0 < 1) or (k0 > self.get_num_ports()):
                raise ValueError("not a valid port number")
            k0se = self.diff_to_se(k0)
            p_ne = self.get_ne_ports(k0)
            for k1 in p_ne:
                k1se = self.diff_to_se(k1)
                if plot_diff is True:
                    Sp_return_k_mag = np.array(self._results.Src(k0se[0] - 1,
                                               k1se[0] - 1, "ComplexMagnitude"))
                    plt.semilogx(np.array(self._results.frequencies()),
                                 20 * np.log10(Sp_return_k_mag),
                                 label='Sd' + str(k0) + 'd' + str(k1))
                    Sp_return_k_mag = np.array(self._results.Src(k0se[1] - 1,
                                               k1se[1] - 1, "ComplexMagnitude"))
                    plt.semilogx(np.array(self._results.frequencies()),
                                 20 * np.log10(Sp_return_k_mag),
                                 label='Sc' + str(k0) + 'c' + str(k1))
                if plot_mixed is True:
                    Sp_return_k_mag = np.array(self._results.Src(k0se[0] - 1,
                                               k1se[1] - 1, "ComplexMagnitude"))
                    plt.semilogx(np.array(self._results.frequencies()),
                                 20 * np.log10(Sp_return_k_mag),
                                 label='Sd' + str(k0) + 'c' + str(k1))
                    Sp_return_k_mag = np.array(self._results.Src(k0se[1] - 1,
                                               k1se[0] - 1, "ComplexMagnitude"))
                    plt.semilogx(np.array(self._results.frequencies()),
                                 20 * np.log10(Sp_return_k_mag),
                                 label='Sc' + str(k0) + 'd' + str(k1))
        plt.legend()
        return fig, ax

    def plot_next_victim(self, port_list, plot_title='Next Victim', plot_diff=True,
                         plot_mixed=False,):
        """Plots the next victim coupling associated with a list of pair of ports. \n
        Depending on the state, 'se' or 'diff', either the single ended
        or diff mode S-parameters are plotted

        Args:
            port_list (list of lists): [[p1, p2], [p3,p4]], all ints
            plot_title (str, optional): Title used for plot. Defaults to 'Next Victim'.
            plot_diff (bool, optional): Used only if state == 'diff' and
            then plots diff to diff and common to common mode components.
            Defaults to True.
            plot_mixed (bool, optional): Used only if state == 'diff' and
            then plots diff to common and common to diff mode components.
            Defaults to False.

        Raises:
            ValueError: Throws an error if an invalid port number is
            specified in port_list

        Returns:
            _type_: fig, ax: Matplotlib objects
        """
        if self.state == "se":
            fig, ax = self.plot_next_victim_se(port_list, plot_title=plot_title)
        elif self.state == "diff":
            fig, ax = self.plot_next_victim_diff(port_list, plot_title=plot_title,
                                                 plot_diff=plot_diff,
                                                 plot_mixed=plot_mixed)
        else:
            raise ValueError("State self.state has an invalid value")
        return fig, ax

    def plot_next_aggressor_se(self, port_list, plot_title='Next Aggressor'):
        """Plots the single ended next aggressor coupling for a list of pair of ports.
        If this function is called and the class is in the 'diff' state,
        the state is changed to 'se'

        Args:
            port_pairs (list of lists): [[p1, p2], [p3,p4]]
            plot_title (str, optional): Title string used in plot.
            Defaults to "Next Aggressor".

        Raises:
            ValueError: Throws an error if an invalid port number is in port_list

        Returns:
           fig, ax: Matplotlib objects
        """
        if self.state == "diff":
            self.compute_se_data()
        fig, ax = plt.subplots(1, 1)
        plt.xlabel('F (Hz)')
        plt.ylabel('dB')
        plt.grid()
        plt.title(plot_title)
        for k0 in port_list:
            p_ne = self.get_ne_ports(k0)
            for k1 in p_ne:
                S_k1_k0 = np.array(self._results.Src(k1 - 1, k0 - 1,
                                   "ComplexMagnitude"))
                plt.semilogx(self._results.frequencies(),
                             20 * np.log10(S_k1_k0),
                             label='S' + str(k1) + ',' + str(k0))
        plt.legend()
        return fig, ax

    def plot_next_aggressor_diff(self, port_list, plot_diff=True, plot_mixed=False,
                                 plot_title='Next Aggressor Diff'):
        """Plots the mixed mode, diff and common, next aggressor coupling
        for a list of pair of ports.
        If this function is called and the class is in the 'se' state,
        the state is changed to 'diff'

        Args:
            diff_port_list (list of lists): [[p1, p2], [p3,p4]]
            plot_diff (True or False): plot the diff and common mode terms.
            Defaults to True
            plot_mixed (True or False): plot the mode conversion terms,
            diff to common and common to diff. Defaults to False
            plot_title (str, optional): Title string used in plot.
            Defaults to "Next Aggressor Diff".

        Raises:
            ValueError: Throws an error if an invalid port number is in port_list

        Returns:
            fig, ax: Matplotlib objects
        """
        if self.state == "se":
            self.compute_differential_data()
        fig, ax = plt.subplots(1, 1)
        plt.xlabel('F (Hz)')
        plt.ylabel('dB')
        plt.grid()
        plt.title(plot_title)
        for k0 in port_list:
            if (k0 < 1) or (k0 > self.get_num_ports()):
                raise ValueError("not a valid port number")
            k0se = self.diff_to_se(k0)
            p_ne = self.get_ne_ports(k0)
            for k1 in p_ne:
                k1se = self.diff_to_se(k1)
                if plot_diff is True:
                    Sp_return_k_mag = np.array(self._results.Src(k1se[0] - 1,
                                               k0se[0] - 1, "ComplexMagnitude"))
                    plt.semilogx(np.array(self._results.frequencies()),
                                 20 * np.log10(Sp_return_k_mag),
                                 label='Sd' + str(k1) + 'd' + str(k0))
                    Sp_return_k_mag = np.array(self._results.Src(k1se[1] - 1,
                                               k0se[1] - 1, "ComplexMagnitude"))
                    plt.semilogx(np.array(self._results.frequencies()),
                                 20 * np.log10(Sp_return_k_mag),
                                 label='Sc' + str(k1) + 'c' + str(k0))
                if plot_mixed is True:
                    Sp_return_k_mag = np.array(self._results.Src(k1se[0] - 1,
                                               k0se[1] - 1, "ComplexMagnitude"))
                    plt.semilogx(np.array(self._results.frequencies()),
                                 20 * np.log10(Sp_return_k_mag),
                                 label='Sd' + str(k1) + 'c' + str(k0))
                    Sp_return_k_mag = np.array(self._results.Src(k1se[1] - 1,
                                               k0se[0] - 1, "ComplexMagnitude"))
                    plt.semilogx(np.array(self._results.frequencies()),
                                 20 * np.log10(Sp_return_k_mag),
                                 label='Sc' + str(k1) + 'd' + str(k0))
        plt.legend()
        return fig, ax

    def plot_next_aggressor(self, port_list, plot_title='Next Aggressor',
                            plot_diff=True, plot_mixed=False,):
        """Plots the next aggressor coupling associated with a list of pair of ports.
        Depending on the state, 'se' or 'diff', either the single ended or
        diff mode S-parameters are plotted

        Args:
            port_list (list of lists): [[p1, p2], [p3,p4]], all ints
            plot_title (str, optional): Title used for plot.
            Defaults to 'Next Aggressor'.
            plot_diff (bool, optional): Used only if state == 'diff' and
            then plots diff to diff and common to common mode components.
            Defaults to True.
            plot_mixed (bool, optional): Used only if state == 'diff' and
            then plots diff to common and common to diff mode components.
            Defaults to False.

        Raises:
            ValueError: Throws an error if an invalid port number is specified in
            port_list

        Returns:
            _type_: fig, ax: Matplotlib objects
        """
        if self.state == "se":
            fig, ax = self.plot_next_aggressor_se(port_list, plot_title=plot_title)
        elif self.state == "diff":
            fig, ax = self.plot_next_aggressor_diff(port_list, plot_title=plot_title,
                                                    plot_diff=plot_diff,
                                                    plot_mixed=plot_mixed)
        else:
            raise ValueError("State self.state has an invalid value")
        return fig, ax

    def plot_fext_victim_se(self, port_list, plot_title='Fext Victim'):
        """Plots the single ended fext victim coupling for a list of pair of ports.
        If this function is called and the class is in the 'diff' state,
        the state is changed to 'se'

        Args:
            port_pairs (list of lists): [[p1, p2], [p3,p4]]
            plot_title (str, optional): Title string used in plot.
            Defaults to "Fext Victim"

        Raises:
            ValueError: Throws an error if an invalid port number is in port_list

        Returns:
           fig, ax: Matplotlib objects
        """
        if self.state == "diff":
            self.compute_se_data()
        fig, ax = plt.subplots(1, 1)
        plt.xlabel('F (Hz)')
        plt.ylabel('dB')
        plt.grid()
        plt.title(plot_title)
        for k0 in port_list:
            if (k0 < 1) or (k0 > self.get_num_ports()):
                raise ValueError("not a valid port number")
            p_ne = self.get_fe_ports(k0)
            for k1 in p_ne:
                S_k0_k1 = np.array(self._results.Src(k0 - 1,
                                   k1 - 1, "ComplexMagnitude"))
                plt.semilogx(self._results.frequencies(),
                             20 * np.log10(S_k0_k1),
                             label='S' + str(k0) + ',' + str(k1))
        plt.legend()
        return fig, ax

    def plot_fext_victim_diff(self, port_list, plot_diff=True, plot_mixed=False,
                              plot_title='Fext Victim Diff'):
        """Plots the mixed mode, diff and common, fext victim coupling
        for a list of pair of ports.
        If this function is called and the class is in the 'se' state,
        the state is changed to 'diff'

        Args:
            diff_port_list (list of lists): [[p1, p2], [p3,p4]]
            plot_diff (True or False): plot the diff and common mode terms.
            Defaults to True
            plot_mixed (True or False): plot the mode conversion terms,
            diff to common and common to diff. Defaults to False
            plot_title (str, optional): Title string used in plot.
            Defaults to "Fext Victim Diff".

        Raises:
            ValueError: Throws an error if an invalid port number is in port_list

        Returns:
            fig, ax: Matplotlib objects
        """
        if self.state == "se":
            self.compute_differential_data()
        fig, ax = plt.subplots(1, 1)
        plt.xlabel('F (Hz)')
        plt.ylabel('dB')
        plt.grid()
        plt.title(plot_title)
        for k0 in port_list:
            if (k0 < 1) or (k0 > self.get_num_ports()):
                raise ValueError("not a valid port number")
            k0se = self.diff_to_se(k0)
            p_ne = self.get_fe_ports(k0)
            for k1 in p_ne:
                k1se = self.diff_to_se(k1)
                if plot_diff is True:
                    Sp_return_k_mag = np.array(self._results.Src(k0se[0] - 1,
                                               k1se[0] - 1, "ComplexMagnitude"))
                    plt.semilogx(np.array(self._results.frequencies()),
                                 20 * np.log10(Sp_return_k_mag),
                                 label='Sd' + str(k0) + 'd' + str(k1))
                    Sp_return_k_mag = np.array(self._results.Src(k0se[1] - 1,
                                               k1se[1] - 1, "ComplexMagnitude"))
                    plt.semilogx(np.array(self._results.frequencies()),
                                 20 * np.log10(Sp_return_k_mag),
                                 label='Sc' + str(k0) + 'c' + str(k1))
                if plot_mixed is True:
                    Sp_return_k_mag = np.array(self._results.Src(k0se[0] - 1,
                                               k1se[1] - 1, "ComplexMagnitude"))
                    plt.semilogx(np.array(self._results.frequencies()),
                                 20 * np.log10(Sp_return_k_mag),
                                 label='Sd' + str(k0) + 'c' + str(k1))
                    Sp_return_k_mag = np.array(self._results.Src(k0se[1] - 1,
                                               k1se[0] - 1, "ComplexMagnitude"))
                    plt.semilogx(np.array(self._results.frequencies()),
                                 20 * np.log10(Sp_return_k_mag),
                                 label='Sc' + str(k0) + 'd' + str(k1))
        plt.legend()
        return fig, ax

    def plot_fext_victim(self, port_list, plot_title='Fext Victim', plot_diff=True,
                         plot_mixed=False,):
        """Plots the next victim coupling associated with a list of pair of ports.
        Depending on the state, 'se' or 'diff', either the
        single ended or diff mode S-parameters are plotted

        Args:
            port_list (list of lists): [[p1, p2], [p3,p4]], all ints
            plot_title (str, optional): Title used for plot.
            Defaults to 'Fext Victim'.
            plot_diff (bool, optional): Used only if state == 'diff' and
            then plots diff to diff and common to common mode components.
            Defaults to True.
            plot_mixed (bool, optional): Used only if state == 'diff' and
            then plots diff to common and common to diff mode components.
            Defaults to False.

        Raises:
            ValueError: Throws an error if an invalid port number is specified
            in port_list

        Returns:
            _type_: fig, ax: Matplotlib objects
        """
        if self.state == "se":
            fig, ax = self.plot_fext_victim_se(port_list, plot_title=plot_title)
        elif self.state == "diff":
            fig, ax = self.plot_fext_victim_diff(port_list, plot_title=plot_title,
                                                 plot_diff=plot_diff,
                                                 plot_mixed=plot_mixed)
        else:
            raise ValueError("State self.state has an invalid value")
        return fig, ax

    def plot_fext_aggressor_se(self, port_list, plot_title='Fext Aggresor'):
        """Plots the single ended fext aggressor coupling for a list of pair of ports.
        If this function is called and the class is in the 'diff' state,
        the state is changed to 'se'

        Args:
            port_pairs (list of lists): [[p1, p2], [p3,p4]]
            plot_title (str, optional): Title string used in plot.
            Defaults to "Fext Aggressor".

        Raises:
            ValueError: Throws an error if an invalid port number is in port_list

        Returns:
           fig, ax: Matplotlib objects
        """
        if self.state == "diff":
            self.compute_se_data()
        fig, ax = plt.subplots(1, 1)
        plt.xlabel('F (Hz)')
        plt.ylabel('dB')
        plt.grid()
        plt.title(plot_title)
        for k0 in port_list:
            if (k0 < 1) or (k0 > self.get_num_ports()):
                raise ValueError("not a valid port number")
            p_ne = self.get_fe_ports(k0)
            for k1 in p_ne:
                S_k0_k1 = np.array(self._results.Src(k1 - 1,
                                   k0 - 1, "ComplexMagnitude"))
                plt.semilogx(self._results.frequencies(), 20 * np.log10(S_k0_k1),
                             label='S' + str(k1) + ',' + str(k0))
        plt.legend()
        return fig, ax

    def plot_fext_aggressor_diff(self, port_list, plot_diff=True, plot_mixed=False,
                                 plot_title='Fext Aggressor Diff'):
        """Plots the mixed mode, diff and common, fext aggressor
        coupling for a list of pair of ports.
        If this function is called and the class is in the 'se' state,
        the state is changed to 'diff'

        Args:
            diff_port_list (list of lists): [[p1, p2], [p3,p4]]
            plot_diff (True or False): plot the diff and common mode terms.
            Defaults to True
            plot_mixed (True or False): plot the mode conversion terms,
            diff to common and common to diff.
            Defaults to False
            plot_title (str, optional): Title string used in plot.
            Defaults to "Fext Aggressor Diff".

        Raises:
            ValueError: Throws an error if an invalid port number is in port_list

        Returns:
            fig, ax: Matplotlib objects
        """
        if self.state == "se":
            self.compute_differential_data()
        fig, ax = plt.subplots(1, 1)
        plt.xlabel('F (Hz)')
        plt.ylabel('dB')
        plt.grid()
        plt.title(plot_title)
        for k0 in port_list:
            if (k0 < 1) or (k0 > self.get_num_ports()):
                raise ValueError("not a valid port number")
            k0se = self.diff_to_se(k0)
            p_ne = self.get_fe_ports(k0)
            for k1 in p_ne:
                k1se = self.diff_to_se(k1)
                if plot_diff is True:
                    Sp_return_k_mag = np.array(self._results.Src(k1se[0] - 1,
                                               k0se[0] - 1, "ComplexMagnitude"))
                    plt.semilogx(np.array(self._results.frequencies()),
                                 20 * np.log10(Sp_return_k_mag),
                                 label='Sd' + str(k1) + 'd' + str(k0))
                    Sp_return_k_mag = np.array(self._results.Src(k1se[1] - 1,
                                               k0se[1] - 1, "ComplexMagnitude"))
                    plt.semilogx(np.array(self._results.frequencies()),
                                 20 * np.log10(Sp_return_k_mag),
                                 label='Sc' + str(k1) + 'c' + str(k0))
                if plot_mixed is True:
                    Sp_return_k_mag = np.array(self._results.Src(k1se[0] - 1,
                                               k0se[1] - 1, "ComplexMagnitude"))
                    plt.semilogx(np.array(self._results.frequencies()),
                                 20 * np.log10(Sp_return_k_mag),
                                 label='Sd' + str(k1) + 'c' + str(k0))
                    Sp_return_k_mag = np.array(self._results.Src(k1se[1] - 1,
                                               k0se[0] - 1, "ComplexMagnitude"))
                    plt.semilogx(np.array(self._results.frequencies()),
                                 20 * np.log10(Sp_return_k_mag),
                                 label='Sc' + str(k1) + 'd' + str(k0))
        plt.legend()
        return fig, ax

    def plot_fext_aggressor(self, port_list, plot_title='Fext Aggressor',
                            plot_diff=True, plot_mixed=False,):
        """Plots the fext aggressor coupling associated with a list of pair of ports.
        Depending on the state, 'se' or 'diff', either the single ended or
        diff mode S-parameters are plotted

        Args:
            port_list (list of lists): [[p1, p2], [p3,p4]], all ints
            plot_title (str, optional): Title used for plot.
            Defaults to 'Fext Aggressor'.
            plot_diff (bool, optional): Used only if state == 'diff' and
            then plots diff to diff and common to common mode components.
            Defaults to True.
            plot_mixed (bool, optional): Used only if state == 'diff' and
            then plots diff to common and common to diff mode components.
            Defaults to False.

        Raises:
            ValueError: Throws an error if an invalid port number is specified in
            port_list

        Returns:
            _type_: fig, ax: Matplotlib objects
        """
        if self.state == "se":
            fig, ax = self.plot_fext_aggressor_se(port_list, plot_title=plot_title)
        elif self.state == "diff":
            fig, ax = self.plot_fext_aggressor_diff(port_list, plot_title=plot_title,
                                                    plot_diff=plot_diff,
                                                    plot_mixed=plot_mixed)
        else:
            raise ValueError("State self.state has an invalid value")
        return fig, ax
