module edo_task_2

using GLMakie
using UUIDs
using LaTeXStrings
using Roots

# Define the solution function
function y_solution(x, a, b, h)
    return a * cosh((x - b) / a) + h
end

import Base.Math.cosh  # Importing the cosh function

"""
    f(x, a, b, h)

Calculate the value of the equation f(x) = cosh((x - b) / a) + h.

# Arguments
- `x`: The function domain values.
- `a`: The amplitude of the function and the scale of the diference between x - b.
- `b`: Horizontal displacement.
- `h`: Vertical displacement.

# Returns
The calculated value of the equation f(x) in x.

"""
f = (x, a, b, h) -> cosh((x - b) / a) + h  


function plot_function(a, b, h)

    fig = Figure(size = (900, 600))
    ax = Axis(
        fig[1, 1],  
        limits = (-10, 10, 0, 1000),
        title = L"f(x)=cosh(\frac{x-b}{a}) + h",
        xlabel = "x",
        ylabel = "y"
    )

    x_values = -9.0:00.1:9.0
    y_values = map( x -> f(x, a, b, h), x_values)
    scatter!(b, f(b, a, b, h), color = :red, markersize = 20, label = "(b, 0)")
    lines!(ax, x_values, y_values, label = "cosh")
    text!(b + 0.5, 20, text = L"(%$a,%$h)", color = :black)
    text!(-4.5, 500, text =  L"f(x)=cosh(\frac{x-b}{a}) + h", color = :black)
    vlines!(ax, [b], color = :red, linestyle = :dash)
    hlines!(ax, [h], xmax = [1], color = :blue, linestyle = :dash)
    display(fig)
    uuid = UUIDs.uuid4()
    save("$uuid.png", fig)
end  


function equations_system(x, ci0, ci1)
    a, b, h = x

    f0 = cosh((ci0[0] - b) / a) + h - ci0[1] # f(0) = 3
    f1 = cosh((ci1[0] - b) / a) + h - ci1[1] # f(1) = 2

    return [f0, f1]
end

function equations_system_solver()
    # Initial guess for the parameters
    x0 = [1.0, 0.0, 0.0]
    equations_system(x)
    # Solve the system of equations
    solution = fzero(equations, x0)
    a, b, h = solution
end


function plot_a_parameter_variation(;
    a_values::Vector{Float64} = fill(1.0, 4),
    b_values::Vector{Float64} = fill(0.0, 4),
    h_values::Vector{Float64} = fill(0.0, 4) )

    # The Grid layout creation, six quadrants
    fig = Figure(size = (1800, 1000))
    quadrant1 = fig[1, 1] = GridLayout()
    quadrant2 = fig[1, 2] = GridLayout()
    quadrant3 = fig[2, 1] = GridLayout()
    quadrant4 = fig[2, 2] = GridLayout()

    a_1 = a_values[1]
    a_2 = a_values[2]
    a_3 = a_values[3]
    a_4 = a_values[4]

    b_1 = b_values[1]
    b_2 = b_values[2]
    b_3 = b_values[3]
    b_4 = b_values[4]

    h_1 = h_values[1]
    h_2 = h_values[2]
    h_3 = h_values[3]
    h_4 = h_values[4]

    quadrants = [
        Axis(
            quadrant1[1, 1], 
            xlabel = L"x", 
            ylabel = L"f(x)",  
            limits = (-10, 10, 0, 1000),
            title = L"f(x)={%$a_1}cosh(\frac{x-%$b_1}{%$a_1}) + %$h_1 ~|~a=%$a_1,~h=%$h_1 ~ y~b=%$b_1"
        ),
        Axis(
            quadrant2[1, 1], 
            xlabel = L"x", 
            ylabel = L"f(x)",  
            limits = (-10, 10, 0, 1000),
            title = L"f(x)={%$a_2}cosh(\frac{x-%$b_2}{%$a_2}) + %$h_2 ~|~a=%$a_2,~h=%$h_2 ~y~b=%$b_2"
        ),
        Axis(
            quadrant3[1, 1], 
            xlabel = L"x", 
            ylabel = L"f(x)",  
            limits = (-10, 10, 0, 1000),
            title = L"f(x)={%$a_3}cosh(\frac{x-%$b_3}{%$a_3}) + %$h_3 ~|~a=%$a_3,~h=%$h_3 ~y~b=%$b_3"
        ),
        Axis(
            quadrant4[1, 1], 
            xlabel = L"x", 
            ylabel = L"f(x)",  
            limits = (-10, 10, 0, 1000),
            title = L"f(x)={%$a_4}cosh(\frac{x-%$b_4}{%$a_4}) + %$h_4 ~|~a=%$a_4,~h=%$h_4 ~y~b=%$b_4"
        )
    ]

    for (index, quadrant) in enumerate(quadrants)
        x_values = -9.0:00.1:9.0
        y_values = map( x -> f(x, a_values[index], b_values[index], h_values[index]), x_values)
        lines!(
            quadrant, 
            x_values, 
            y_values, 
            label = "cosh")
        a_tag = a_values[index]
        h_tag = h_values[index]
        vlines!(
            quadrant, 
            [b_values[index]], 
            color = :red, 
            linestyle = :dash)
        hlines!(
            quadrant, 
            [h_values[index]], 
            xmax = [1], 
            color = :blue, 
            linestyle = :dash)          
    end
    uuid = UUIDs.uuid4()
    save("$uuid.png", fig)  
    display(fig)
end    


function plot_electric_field(a, b, q)
    f(x, y) = (q/2) * log((x-a)^2 + (y-b)^2)
    xs = range(-15, stop = 15, length = 1000)
    ys = range(-15, stop = 15, length = 1000)
    #z = f.(x', y)
    zs = [f(x, y) for x in xs, y in ys]
   
end

function gradient_p(x, y, a, b, q)
    #z = f.(x', y)
    denominator = ((x - a)^2) + ((y - b)^2)
    dp_dx = q * (x - a) / denominator
    dp_dy = q * (y - b) / denominator
    return Point2(dp_dx, dp_dy) 
end

function gradient_p_w(a, b, q)
    xs = range(-5, stop = 5, length = 20)
    ys = range(-5, stop = 5, length = 20)
    zs = [gradient_p(x, y, a, b, q) for x in xs, y in ys]
    z =  reduce(vcat, collect(eachrow(zs)))
    
    println(z)
    return z
end

function dp_dx(x, y, a, b, q)
    return q * (x - a) / ((x - a)^2 + (y - b)^2)
end

function dp_dy(x, y, a, b, q)
    return q * (y - b) / ((x - a)^2 + (y - b)^2)
end

function inverse_p(a, b, q)
    f(x, y) = a + sqrt( abs(exp(2x / q) - (y - b)^2))
    xs = range(-10, stop = 10, length = 1000)
    ys = range(-10, stop = 10, length = 1000)
    #z = f.(x', y)
    zs = [f(x, y) for x in xs, y in ys]
    return zs
end

function electric_field()

    fig = Figure()
    ax1 = Axis(fig[1, 1])
    qs = [(0.0, 1.0, -1.0), (0.0, -1.0, -1.0), (1.0, 0.0, -1.0)]
    #interations = map(q -> plot_electric_field(q[1], q[2], q[3]), qs) 
    xs = range(-2.5, stop = 2.5, length = 2)
    ys = range(-2.5, stop = 2.5, length = 2)

    
    #contour!( foldl(+, interations), levels=-1:0.1:1, labels=true )
   
    function es(x, y)
        return reduce( (state, t) -> Point2f( state[1] + t[1], state[2] + t[2]), 
                                    map( q -> gradient_p(x, y, q[1], q[2], q[3]), qs), init = Point2f(0, 0))
    end   
    
    zs = [es(x, y) for x in xs, y in ys]
println(zs)
    
    arrows!( ax1, xs, ys, x -> es(x[1], x[2]), levels=-2:0.1:2, arrowsize = 7, lengthscale = 0.25 )     


    #fels, sels = zip(interactions...)
    

   #xs = range(-15, stop = 15, length = 1000)
    #ys = range(-15, stop = 15, length = 1000)
    #U = map(q -> dp_dx(xs, ys, q[1], q[2],[3]), qs)
    #V = map(q -> dp_dy(xs, ys, q[1], q[2], q[3]), qs)
    #fig, ax, pl = streamplot(ff, -1.5..1.5, -1.5..1.5, colormap = :magma)
    #streamplot!(ax2, ff, -1.5 .. 1.5, -1.5 .. 1.5, color = :blue)
    #for q in qs
        #streamplot!(fig[1,1], (x, y) -> gradient_p(x, y, q[1], q[2], q[3]), -1.5 .. 1.5, -1.5 .. 1.5, colormap = :magma)
    #end
    display(fig)
end


function  alfo()
    

    # Set parameters
    a = 1.0
    b = 0.0
    h = 1.0

    # Generate x values
    x_values = range(-5, stop=5, length=100)

    # Calculate y values using the solution function
    y_values = y_solution.(x_values, a, b, h)

    # Create Makie plot
    fig = Figure()
    ax = Axis(fig[1, 1])
    lines!(ax, x_values, y_values, label="y(x)")
    save("output_plot.png", fig)
    display(fig)
end

function pendulum_solver1(theta0, theta_dot0, g, l, dt, num_steps)
    # Initialize arrays to store the values of theta and theta_dot
    theta = [theta0]
    theta_dot = [theta_dot0]

    # Perform the numerical integration
    for i in 1:num_steps
        # Compute the four intermediate slopes (k1, k2, k3, k4)
        k1 = dt * theta_dot[i]
        k2 = dt * (theta_dot[i] + k1/2)
        k3 = dt * (theta_dot[i] + k2/2)
        k4 = dt * (theta_dot[i] + k3)

        # Update theta and theta_dot using the fourth-order Runge-Kutta method
        theta_new = theta[i] + (k1 + 2*k2 + 2*k3 + k4) / 6
        theta_dot_new = theta_dot[i] - (g / l) * sin(theta[i]) * dt

        # Append the new values to the arrays
        push!(theta, theta_new)
        push!(theta_dot, theta_dot_new)
    end

    return theta, theta_dot
end

# Set the initial conditions and parameters
#= theta0 = 0.1  # Initial angular displacement
theta_dot0 = 0.0  # Initial angular velocity
g = 9.81  # Acceleration due to gravity
l = 1.0  # Length of the pendulum
dt = 0.01  # Time step size
num_steps = 1000  # Number of iterations

# Solve the nonlinear pendulum differential equation
theta, theta_dot = pendulum_solver1(theta0, theta_dot0, g, l, dt, num_steps)

# Print the values of theta at each time step
for i in 1:num_steps
    println("Theta[$i]: ", theta[i])
end

# Plot the trajectory of the pendulum
plot(theta, label="Theta")

# Save the plot
save("pendulum_trajectory.png") =#

# Display the plot
# display(plot!())

function pendulum_solver2(theta0, theta_dot0, g, l, dt, num_steps)
    # Initialize arrays to store the values of theta and theta_dot
    theta = [theta0]
    theta_dot = [theta_dot0]

    # Perform the numerical integration
    for i in 1:num_steps
        # Compute the four intermediate slopes (k1, k2, k3, k4)
        k1 = dt * theta_dot[i]
        k2 = dt * (theta_dot[i] + k1/2)
        k3 = dt * (theta_dot[i] + k2/2)
        k4 = dt * (theta_dot[i] + k3)

        # Update theta and theta_dot using the fourth-order Runge-Kutta method
        theta_new = theta[i] + (k1 + 2*k2 + 2*k3 + k4) / 6
        theta_dot_new = theta_dot[i] - (g / l) * sin(theta[i]) * dt

        # Append the new values to the arrays
        push!(theta, theta_new)
        push!(theta_dot, theta_dot_new)
    end

    return theta, theta_dot
end

# Set the parameters
#= g = 9.81  # Acceleration due to gravity
l = 1.0  # Length of the pendulum
dt = 0.01  # Time step size
num_steps = 1000  # Number of iterations

# Define the initial conditions for three different theta values
initial_conditions = [0.1, 0.5, 1.0]

# Loop over the initial conditions and solve the nonlinear pendulum differential equation
for (idx, theta0) in enumerate(initial_conditions)
    # Set the initial angular velocity to zero for all cases
    theta_dot0 = 0.0

    # Solve the nonlinear pendulum differential equation
    theta, theta_dot = pendulum_solver2(theta0, theta_dot0, g, l, dt, num_steps)

    # Print the values of theta at each time step
    println("Results for initial condition theta0 = $theta0:")
    for i in 1:num_steps
        println("Theta[$i]: ", theta[i])
    end
    println()
end =#


function pendulum_solver3(theta0, theta_dot0, g, l, dt, num_steps)
    # Initialize arrays to store the values of theta and theta_dot
    theta = [theta0]
    theta_dot = [theta_dot0]

    # Perform the numerical integration
    for i in 1:num_steps
        # Compute the four intermediate slopes (k1, k2, k3, k4)
        k1 = dt * theta_dot[i]
        k2 = dt * (theta_dot[i] + k1/2)
        k3 = dt * (theta_dot[i] + k2/2)
        k4 = dt * (theta_dot[i] + k3)

        # Update theta and theta_dot using the fourth-order Runge-Kutta method
        theta_new = theta[i] + (k1 + 2*k2 + 2*k3 + k4) / 6
        theta_dot_new = theta_dot[i] - (g / l) * sin(theta[i]) * dt

        # Append the new values to the arrays
        push!(theta, theta_new)
        push!(theta_dot, theta_dot_new)
    end

    return theta, theta_dot
end

# Set the parameters


function plottingNonNumerialAproximation()
    g = 9.81  # Acceleration due to gravity
    l = 1.0  # Length of the pendulum
    dt = 0.01  # Time step size
    num_steps = 1000  # Number of iterations

    # Set the initial conditions
    theta0 = 0.1  # Initial angular displacement
    theta_dot0 = 0.0  # Initial angular velocity

    x_values = range(1, stop=1000, length=1001)
    # Solve the nonlinear pendulum differential equation
    theta, _ = pendulum_solver3(theta0, theta_dot0, g, l, dt, num_steps)

    # Plot the graph using GLMakie
   #=  scene = Scene()
    lines!(scene, theta, color=:blue, linewidth=1)
    #title!(scene, "Nonlinear Pendulum Motion")
    display(scene) =#
    fig = Figure()
    ax = Axis(fig[1, 1])
    lines!(ax, x_values, theta, label="y(x)")
    uuid = UUIDs.uuid4()
    save("$uuid.png", fig)
    display(fig)
end

end # module edo_task_2
