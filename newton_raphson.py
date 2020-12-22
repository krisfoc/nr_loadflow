import numpy as np


class Bus:
    def __init__(self, Busdata):
        self.num = int(Busdata[0])
        self.type = Busdata[1]
        self.P_spec = Busdata[2]
        self.Q_spec = Busdata[3]
        self.V_spec = Busdata[4]
        self.Qlim = [Busdata[5], Busdata[6]]
        self.Plim = [Busdata[7], Busdata[8]]
        if not np.isnan(Busdata[9]):
            self.comp = Busdata[9]
        else:
            self.comp = 0

        self.connected_lines = {}  # type: Dict[Bus,Line]
        self.powerflow = {}  # type: Dict[Bus,float]
        # Load flow internal and result variables
        self.n = 0
        self.angle = 0
        if not np.isnan(self.V_spec):
            self.V = self.V_spec
        else:
            self.V = 1
        self.P = 0
        self.Q = 0
        if Busdata[1] == 'PV':
            self.is_PV = True
        else:
            self.is_PV = False

    def print(self):

        stri = "{:^6}|{:^6}|{:^8}|{:^8}|{:^8}|{:^8}|".format(self.num, self.type, self.P_spec, self.Q_spec, self.V_spec,
                                                             np.round(self.angle, 4))
        stri += "{:^8}|{:^8}|{:^8}".format(np.round(self.V, 4), np.round(self.P, 4), np.round(self.Q, 4))
        print(stri)

    def printpf(self):
        for bus, pf in self.powerflow.items():
            print('From {} to {}:  {}'.format(self.num, bus.num, np.round(pf, 4)))


class Line:
    def __init__(self, Linedata):
        self.nodes = Linedata[0:2]
        self.R = Linedata[2]
        self.X = Linedata[3]
        self.y = 1 / (self.R + 1j * self.X)

    def other(self, bus):
        if bus.num != self.nodes[1]:
            return self.nodes[1]
        else:
            return self.nodes[0]


def Qc(i, j, y):
    # returns -(Gij*sin(ti-tj)-Bij*cos(ti-tj))
    return y.real * np.sin(i.angle - j.angle) - y.imag * np.cos(i.angle - j.angle)


def Qcdj(i, j, y):
    # returns -(Gij*sin(ti-tj)-Bij*cos(ti-tj))/dtj
    return -y.real * np.cos(i.angle - j.angle) - y.imag * np.sin(i.angle - j.angle)


def Qcdi(i, j, y):
    # returns -(Gij*sin(ti-tj)-Bij*cos(ti-tj))/dti
    return y.real * np.cos(i.angle - j.angle) + y.imag * np.sin(i.angle - j.angle)


def Pc(i, j, y):
    # returns -(Gij*cos(ti-tj)+Bij*sin(ti-tj))
    return y.real * np.cos(i.angle - j.angle) + y.imag * np.sin(i.angle - j.angle)


def Pcdj(i, j, y):
    # -(Gij*cos(ti-tj)+Bij*sin(ti-tj))/dtj
    return y.real * np.sin(i.angle - j.angle) - y.imag * np.cos(i.angle - j.angle)


def Pcdi(i, j, y):
    # -(Gij*cos(ti-tj)+Bij*sin(ti-tj))/dti
    return -y.real * np.sin(i.angle - j.angle) + y.imag * np.cos(i.angle - j.angle)


def calculate_powers(bus_lst):
    for bus1 in bus_lst:
        bus1.P = 0
        bus1.Q = 0
        total_y = 1j * bus1.comp
        for bus2 in bus1.connected_lines:
            total_y += bus1.connected_lines[bus2].y
            bus1.P -= bus1.V * bus2.V * Pc(bus1, bus2, bus1.connected_lines[bus2].y)
            bus1.Q -= bus1.V * bus2.V * Qc(bus1, bus2, bus1.connected_lines[bus2].y)
        bus1.P += bus1.V ** 2 * total_y.real
        bus1.Q -= bus1.V ** 2 * total_y.imag

    return

#prints all the data of the system
def printsys(bus_lst):
    str = "{:^6}|{:^6}|{:^8}|{:^8}|{:^8}|{:^8}|{:^8}|{:^8}|{:^8}".format("#", "Type", "Pspec", "Qspec",
                                                                         "Vspec", "theta", "V", "P", "Q")
    print(str)
    for bus in bus_lst:
        bus.print()
    return

#returns the numerator of the other bus on the same line
def get_other_line(i, line):
    if i == line.nodes[0]:
        return line.nodes[1]
    else:
        return line.nodes[0]


#returns the relevant value the of entry in the jacobian
def derivate(bus_lst, li1, li2):
    val = 0
    this = li1[1]
    if li1[0] == 'p':
        if li2[0] == 'v':
            if li1[1] == li2[1]:
                val = 0
                g = 0
                for item, line in bus_lst[this].connected_lines.items():
                    other = get_other_line(this, line)
                    g += line.y.real
                    val -= bus_lst[other].V * Pc(bus_lst[this], bus_lst[other], line.y)
                val += 2 * bus_lst[this].V * g
                return val
            else:
                val = 0
                for item, line in bus_lst[this].connected_lines.items():
                    other = get_other_line(this, line)
                    if other == li2[1]:
                        val -= bus_lst[this].V * Pc(bus_lst[this], bus_lst[other], line.y)
                return val
        if li2[0] == 't':
            if li1[1] == li2[1]:
                val = 0
                for item, line in bus_lst[this].connected_lines.items():
                    other = get_other_line(this, line)
                    val -= bus_lst[this].V * bus_lst[other].V * Pcdi(bus_lst[this], bus_lst[other], line.y)
                return val
            else:
                val = 0
                for item, line in bus_lst[this].connected_lines.items():
                    other = get_other_line(this, line)
                    if other == li2[1]:
                        val -= bus_lst[this].V * bus_lst[other].V * Pcdj(bus_lst[this], bus_lst[other], line.y)
                return val

    if li1[0] == 'q':
        if li2[0] == 'v':
            if li1[1] == li2[1]:
                val = 0
                b = bus_lst[this].comp
                for item, line in bus_lst[this].connected_lines.items():
                    other = get_other_line(this, line)
                    b += line.y.imag
                    val -= bus_lst[other].V * Qc(bus_lst[this], bus_lst[other], line.y)
                val -= 2 * bus_lst[this].V * b
                return val
            else:
                val = 0
                for item, line in bus_lst[this].connected_lines.items():
                    other = get_other_line(this, line)
                    if other == li2[1]:
                        val -= bus_lst[this].V * Qc(bus_lst[this], bus_lst[other], line.y)
                return val
        if li2[0] == 't':
            if li1[1] == li2[1]:
                val = 0
                for item, line in bus_lst[this].connected_lines.items():
                    other = get_other_line(this, line)
                    val -= bus_lst[this].V * bus_lst[other].V * Qcdi(bus_lst[this], bus_lst[other], line.y)
                return val
            else:
                val = 0
                for item, line in bus_lst[this].connected_lines.items():
                    other = get_other_line(this, line)
                    if other == li2[1]:
                        val -= bus_lst[this].V * bus_lst[other].V * Qcdj(bus_lst[this], bus_lst[other], line.y)
                return val

    return val


def calculate_jacobian(bus_lst):
    jac = []
    # x = known powers
    x = []
    # w = state variables
    w = []
    t_er = []
    v_er = []
    p_er = []
    q_er = []
    for bus in bus_lst:
        if bus.type == 'PQ':
            p_er.append(['p', bus.num])
            q_er.append(['q', bus.num])
            v_er.append(['v', bus.num])
            t_er.append(['t', bus.num])
        elif bus.type == 'PV':
            p_er.append(['p', bus.num])
            t_er.append(['t', bus.num])
        else:
            None
    x.extend(p_er)
    x.extend(q_er)
    w.extend(t_er)
    w.extend(v_er)
    # print(x)
    # print(w)
    length = len(x)
    for i in range(0, length):
        temp = []
        for j in range(0, length):
            temp.append(derivate(bus_lst, x[i], w[j]))
        jac.append(temp)
    # print(jac)
    return jac, x, w


def newton_raphson_iteration(bus_lst, jac, x, w):
    dP = []
    dQ = []
    dx = []
    for item in x:
        for bus in bus_lst:
            if bus.num == item[1] and item[0] == 'p':
                dP.append(bus.P_spec - bus.P)
            if bus.num == item[1] and item[0] == 'q':
                dQ.append(bus.Q_spec - bus.Q)
    dx.extend(dP)
    dx.extend(dQ)
    dw = np.linalg.solve(jac, dx)
    # print('w:', w)
    # print('dw:', dw)
    for c in range(0, len(w)):
        for bus in bus_lst:
            if w[c][1] == bus.num:
                if w[c][0] == 'v':
                    bus.V += dw[c]
                else:
                    bus.angle += dw[c]
    return


# need to make check to switch previous PV-buses back to PV
def pv_limits(bus_lst):
    for bus in bus_lst:
        if bus.type == 'PV':
            if bus.Q > bus.Qlim[1]:
                print('Upper Q limit violated at bus ', bus.num)
                bus.type = 'PQ'
                bus.Q_spec = bus.Qlim[1]
            if bus.Q < bus.Qlim[0]:
                print('Lower Q limit violated at bus ', bus.num)
                bus.type = 'PQ'
                bus.Q_spec = bus.Qlim[0]
        if bus.type == 'PQ' and bus.is_PV:
            if bus.V > bus.V_spec:
                None
    return


# returns the voltage as a phasor
def phasor(bus):
    return bus.V * np.exp(bus.angle * 1j)


def powerflows(bus_lst):
    injected_power = 0
    recieved_power = 0
    print('transmission line power flow')
    print(' Bus # to #       P  +  Q    ')
    for bus1 in bus_lst:
        for bus2, line in bus1.connected_lines.items():
            I = (phasor(bus1) - phasor(bus2)) * line.y
            bus1.powerflow[bus2] = phasor(bus1) * np.conj(I)
            if bus1.powerflow[bus2] >= 0:
                injected_power += bus1.powerflow[bus2]
            else:
                recieved_power += bus1.powerflow[bus2]
        bus1.printpf()
    losses = injected_power + recieved_power
    print('power losses:', np.around(losses, 5))


def newton_raphson(bus_lst, set_print=True):
    calculate_powers(bus_lst)

    if set_print and False: #this if does not run
        print('---------------------------------------------------------------------------')
        print('initial conditions')
        printsys(bus_lst)
        # print('Initial Jacobian')
        jac, x, w = calculate_jacobian(bus_lst)
        # print(jac)
    # should converge after 4 or 5 iteration
    for i in range(4):
        pv_limits(bus_lst)
        jac, x, w = calculate_jacobian(bus_lst)
        # print(jac)
        newton_raphson_iteration(bus_lst, jac, x, w)
        calculate_powers(bus_lst)

        if set_print:
            print('---------------------------------------------------------------------------')
            print('iteration:', i + 1)
            printsys(bus_lst)
    if set_print:
        powerflows(bus_lst)
        None

    return
