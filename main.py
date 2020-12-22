import newton_raphson as nr

#example
def newton_raphson_2():
    z12 = 0.05 + 0.25j
    z23 = 0.05 + 0.25j
    z23 = z23 * z23 / (z23 + z23)
    z14 = 0.02 + 0.2j
    z24 = 0.02 + 0.2j
    z45 = 0.01 + 0.1j
    p1 = 1
    p2 = -0.6
    q2 = -0.3
    p4 = -0.6
    q4 = -0.2
    p5 = -0.5
    q5 = -0.4
    # [num, type, P_spec, Q_spec, V_spec, Qlim, Qlim, Plim, Plim, comp]
    bus0 = nr.Bus([0, 'PV', p1, 0, 1, -100, 100, -100, 100, 0])
    bus1 = nr.Bus([1, 'PQ', p2, q2, 1, -100, 100, -100, 100, 0])
    bus2 = nr.Bus([2, 'SB', 0, 0, 1, -100, 100, -100, 100, 0])
    bus3 = nr.Bus([3, 'PQ', p4, q4, 1, -100, 100, -100, 100, 0])
    bus4 = nr.Bus([4, 'PQ', p5, q5, 1, -100, 100, -100, 100, 0])

    # [from, to, R, X]
    Line1 = nr.Line([0, 1, z12.real, z12.imag])
    Line2 = nr.Line([0, 3, z14.real, z14.imag])
    Line3 = nr.Line([1, 3, z24.real, z24.imag])
    Line4 = nr.Line([1, 2, z23.real, z23.imag])
    Line5 = nr.Line([3, 4, z45.real, z45.imag])

    bus0.connected_lines[bus1] = Line1
    bus1.connected_lines[bus0] = Line1
    bus0.connected_lines[bus3] = Line2
    bus3.connected_lines[bus0] = Line2
    bus1.connected_lines[bus3] = Line3
    bus3.connected_lines[bus1] = Line3
    bus1.connected_lines[bus2] = Line4
    bus2.connected_lines[bus1] = Line4
    bus3.connected_lines[bus4] = Line5
    bus4.connected_lines[bus3] = Line5


    bus_data = [bus0, bus1, bus2, bus3, bus4]
    nr.newton_raphson(bus_data, set_print=True)
    angles = [i.angle for i in bus_data]
    Vs = [i.V for i in bus_data]
    Vs_real_values=[0, 0, 0, 0, 0]
    for i in range(0, 5):
        Vs_real_values[i]=int(132*1000*Vs[i])
    print('real voltage magnitudes: ', Vs_real_values)
    return bus_data


newton_raphson_2()
