from __future__ import division
from itertools import count
from models.Buses import Buses
from scripts.stamp_helpers import *
from models.global_vars import global_vars

class Loads:
    _ids = count(0)

    def __init__(self,
                 Bus,
                 P,
                 Q,
                 IP,
                 IQ,
                 ZP,
                 ZQ,
                 area,
                 status):
        """Initialize an instance of a PQ or ZIP load in the power grid.

        Args:
            Bus (int): the bus where the load is located
            self.P (float): the active power of a constant power (PQ) load.
            Q (float): the reactive power of a constant power (PQ) load.
            IP (float): the active power component of a constant current load.
            IQ (float): the reactive power component of a constant current load.
            ZP (float): the active power component of a constant admittance load.
            ZQ (float): the reactive power component of a constant admittance load.
            area (int): location where the load is assigned to.
            status (bool): indicates if the load is in-service or out-of-service.
        """
        self.Bus = Bus
        self.P_MW = P
        self.Q_MVA = Q
        self.IP_MW = IP
        self.IQ_MVA = IQ
        self.ZP_MW = ZP
        self.ZQ_MVA = ZQ
        self.area = area
        self.status = status
        self.id = Loads._ids.__next__()

        # Convert from MVA to pu
        self.P = P/global_vars.base_MVA
        self.Q = Q/global_vars.base_MVA
        self.IP = IP/global_vars.base_MVA
        self.IQ = IQ/global_vars.base_MVA
        self.ZP = ZP/global_vars.base_MVA
        self.ZQ = ZQ/global_vars.base_MVA

    def assign_indexes(self, bus):
        # Nodes shared by generators on the same bus
        self.Vr_node = bus[Buses.bus_key_[self.Bus]].node_Vr
        self.Vi_node = bus[Buses.bus_key_[self.Bus]].node_Vi

    def stamp(self, V, Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J):
        # self.Vr_node is index within bus voltage array of the ECF real
        # node for this load element, Vi_node is analagous
        Vr = V[self.Vr_node]
        Vi = V[self.Vi_node]

        # Review: Loads are a VCCS, ICS, and Conductance in Parallel
        # ICS: Dependent on previous iteration Irl_hist (real) or Iil_hist (imag)
        # VCCS: Coupled with imaginary dual circuits, eg.
        #   Real VCCS depends on imaginary voltage dIrldVi*Vil and
        #   Imaginary VCSS depends on real voltage dIildVr*Vrl
        # Conductance: dIrldVr*Vr and dIildVi*Vi

        # Linearized Real Equivalent Load Stamps
        Irl_hist = (self.P*Vr+self.Q*Vi)/(Vr**2+Vi**2) # Previous iteration Ir
        dIrldVr = (self.P*(Vi**2-Vr**2) - 2*self.Q*Vr*Vi)/(Vr**2+Vi**2)**2
        dIrldVi = (self.Q*(Vr**2-Vi**2) - 2*self.P*Vr*Vi)/(Vr**2+Vi**2)**2
        Vr_J_stamp = -Irl_hist + dIrldVr*Vr + dIrldVi*Vi

        # Conductance
        idx_Y = stampY(self.Vr_node, self.Vr_node, dIrldVr, Y_val, Y_row, Y_col, idx_Y)

        # VCCS
        idx_Y = stampY(self.Vr_node, self.Vi_node, dIrldVi, Y_val, Y_row, Y_col, idx_Y)

        # ICS
        idx_J = stampJ(self.Vr_node, Vr_J_stamp, J_val, J_row, idx_J)

        # Linearized Real Equivalent Load Stamps
        Iil_hist = (self.P*Vi-self.Q*Vr)/(Vr**2+Vi**2) # Previous iteration Ii
        dIildVi = -dIrldVr
        dIildVr = dIrldVi
        Vi_J_stamp = -Iil_hist + dIildVr*Vr + dIildVi*Vi

        # Conductance
        idx_Y = stampY(self.Vi_node, self.Vr_node, dIildVr, Y_val, Y_row, Y_col, idx_Y)

        # VCCS
        idx_Y = stampY(self.Vi_node, self.Vi_node, dIildVi, Y_val, Y_row, Y_col, idx_Y)

        # ICS
        idx_J = stampJ(self.Vi_node, Vi_J_stamp, J_val, J_row, idx_J)

        return (idx_Y, idx_J)

    def calc_residuals(self, resid, V):
        P = self.P
        Vr = V[self.Vr_node]
        Vi = V[self.Vi_node]
        Q = self.Q
        resid[self.Vr_node] += (P*Vr+Q*Vi)/(Vr**2+Vi**2)
        resid[self.Vi_node] += (P*Vi-Q*Vr)/(Vr**2+Vi**2)
