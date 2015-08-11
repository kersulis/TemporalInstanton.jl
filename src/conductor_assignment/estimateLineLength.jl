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
  lengths vector of line lengths [miles]
"""
function estimateLineLength(System)
# Set default length
default = 1e-4

# Initialize output
lengths = zeros(length(System.Line),1)

for i = 1:length(System.Line)
    if isempty(System.Line[i].conductorModel.R_unit)
        # Line has no conductor model, so set it to the default length
        lengths[i] = default
    else
        # Line has a conductor model, so estimate its length.
        # Find pu base for resistance for this line
        rBase = System.Bus([System.Bus[:].num]==System.Line[i].from).vBase^2/System.Params.sBase

        # Find resistance of line in ohms
        rOhms = System.Line[i].R*rBase

        # Estimate line length in miles (mi = ohms/[ohms/m]/[m/mi])
        lengths[i] = rOhms/(System.Line[i].conductorModel.R_unit/...
            System.Line[i].conductorModel.bundle)/1609.34
    end
end

# Ensure that all length estimates are positive
lengths[lengths.<=0] = default

return lengths
end
