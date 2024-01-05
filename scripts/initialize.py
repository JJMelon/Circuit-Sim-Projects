import numpy as np
def initialize(size_Y, bus, generator, slack, flat_start):
    v_init = np.zeros(size_Y, dtype=np.float) # create zero vector of size_Y

    # For flat start:
    # bus real node voltage is 1 pu
    # bus imaginary node voltage is 0 pu
    # gen reactive power is at average of limits
    # slack current at 0
    if flat_start:
        for ele in bus:
            v_init[ele.node_Vr] = 1
            v_init[ele.node_Vi] = 0
        for ele in generator:
            v_init[ele.Q_node] += (ele.Qmax + ele.Qmin)/2 # Should be just =?
    else:
        for ele in bus:
            v_init[ele.node_Vr] = ele.Vr_init
            v_init[ele.node_Vi] = ele.Vi_init
        for ele in generator:
            v_init[ele.Q_node] += -ele.Qinit
        for ele in slack:
            v_init[ele.Slack_Ir_node] = ele.Ir_init
            v_init[ele.Slack_Ii_node] = ele.Ii_init
    return v_init


