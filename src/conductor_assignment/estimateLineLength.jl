""" 2015-08-11. Ported from .m file provided by Jonathon Martin
This function uses line resistance values and unit resistance values to
estimate line length.
Input
  System  structure containing system information, specifically...
   .Bus
       .num        bus number
       .vBase      base voltage [kV]
   .Line   nline x 1 vector
       .from       from bus number
       .to         to bus number
       .R          resistance [pu]
       .isXfrmr    logical specifying if line is a transformer
       .conductorModel
           .R_unit     unit resistance of conductor [ohms/meter]
   .Params
       .sBase      3-phase power system base [MVA/pu]

Output
  `lengths` vector of line lengths [miles]
"""
function estimateLineLength(ps::PowerSystem)
    # Set default length
    default = 1e-4

    # Initialize output
    lengths = zeros(length(ps.line))

    for i = 1:length(ps.line)
        if isnan(ps.line[i].conductorModel.R_unit)
            # Line has no conductor model, so set it to the default length
            lengths[i] = default
        else
            # Line has a conductor model, so estimate its length.
            # Find pu base for resistance for this line
            rBase = ps.bus[findfirst([b.num for b in ps.bus].==ps.line[i].from)].vBase^2/ps.params.sBase

            # Find resistance of line in ohms
            rOhms = ps.line[i].R*rBase

            # Estimate line length in miles (mi = ohms/[ohms/m]/[m/mi])
            lengths[i] = rOhms/(ps.line[i].conductorModel.R_unit/
                ps.line[i].conductorModel.bundle)# /1609.34
        end
    end

    # Ensure that all length estimates are positive
    lengths[lengths.<=0] = default

    return lengths
end
