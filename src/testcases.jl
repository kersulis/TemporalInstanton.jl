"Test case used for timing analysis"
function testcase(name::ASCIIString)
    if name == "timing"
        d = load_rts96_data(return_as_type=true);
        # Thermal model parameters:
        d.Tamb = 35. # C
        d.T0 = 60. #46. # initial line steady-state temp

        d.time_values = 0:30:300 # five minutes in 30-sec steps
        d.int_length = 300. # seconds = 5 min
        Gp,Dp,Rp = d.G0,d.D0,d.R0
        d.G0 = [0.7*Gp;0.7*Gp;0.7*Gp;0.7*Gp;0.7*Gp;0.7*Gp]
        d.D0 = [0.9*Dp;0.9*Dp;0.9*Dp;0.9*Dp;0.9*Dp;0.9*Dp]
        d.R0 = [Rp;1.1*Rp;1.2*Rp;1.3*Rp;1.4*Rp;1.5*Rp]
        return d
    else
        error("No match.")
    end
end
