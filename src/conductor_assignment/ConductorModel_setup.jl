""" 2015-08-11. Ported from .m file provided by Jonathon Martin
This function uses power rating information to estimate an appropriate
ACSR conductor model to use on each line in the system.

Input
System: structure containing network information

Output
model: structure containing a conductorModel structure for each network
connection.
Fields in conductorModel are
  .R_unit     unit resistance of conductor [ohms/m]
  .I_max      current limit of network connection [A/phase]
  .D          diameter of complete conductor [mm]
  .Al_m       unit mass of aluminum portion of ACSR [kg/m]
  .St_m       unit mass of steel portion of ACSR [kg/m]
  .bundle     number of conductors used in the network connection
"""
function ConductorModel_setup(System)
    for i = 1:System.Params.numOfLine
        if !System.Line[i].isXfrmr && System.Line[i].R > 0
            # determine the voltage base of the line
            vBase = System.Bus[[System.Bus.num] .== System.Line[i].from].vBase

            # determine current limit of the line
            I_max = System.Line[i].sLim*1e3)/(vBase*sqrt(3.)
            if I_max == 0
                # If the calculated rating is zero, assume line is large enough
                # that overloading is not going to happen.
                R_unit = []
                I_max = []
                D = []
                Al_m = []
                St_m = []
                bundle = []
            else
                # match the line to an appropriate ACSR conductor
                D,Al_m,St_m,R_unit,bundle = ACSR_assignment(I_max,vBase)

                if bundle > 6
                    # This network connection most likely has unrealistic flow
                    # limits and should not have a conductor or thermal model.
                    R_unit = []
                    I_max = []
                    D = []
                    Al_m = []
                    St_m = []
                    bundle = []
                end
            end

            # store the results for output
            model.conductorModel[i] = struct('R_unit',R_unit,'I_max',I_max,'D',D,'Al_m',Al_m,'St_m',St_m,'bundle',bundle)
        else
            # if connection is a transformer, store an empty model
            model.conductorModel[i] = struct('R_unit',[],'I_max',[],'D',[],'Al_m',[],'St_m',[],'bundle',[])
        end
    end
    return model
end
