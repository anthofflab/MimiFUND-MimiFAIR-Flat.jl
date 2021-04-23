using Mimi

# import the macro that allows a parameter or variable to access a missing value 
# without throwing an error, which is less safe for users but avoids some corner
# case problems in the context of this type of connecting, unit-conversion component
import Mimi: @allow_missing 

@defcomp multiplier begin

    multiplier  = Parameter()
    input       = Parameter(index=[time])
    output      = Variable(index=[time])

    function run_timestep(p, v, d, t)
        v.output[t] = @allow_missing(p.input[t]) * p.multiplier
    end
end
