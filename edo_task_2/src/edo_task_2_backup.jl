module edo_task_2

# The imported libraries for using in work deveploment.
using GLMakie                   # Plotting framework
using UUIDs                     # UUID Generator
using LaTeXStrings              # For render latex text into the makie plots.
using Roots                     # For linear equation system solving.
using DifferentialEquations     # Numerical resolution for differential equations. 

# Define the solution function
function y_solution(x, a, b, h)
    return a * cosh((x - b) / a) + h    # The catenary function definition.
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
f = (x, a, b, h) -> a*cosh((x - b) / a) + h  


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
    text!(b + 0.5, 20, text = L"(%$b,%$h)", color = :black)
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



function equations_system_solver(x)
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
    qs = [(0.0, 0.5, -1.0), (0.0, -0.5, -1.0), (0.66, 0.0, -1.0)]
    #interations = map(q -> plot_electric_field(q[1], q[2], q[3]), qs) 
    xs = range(-2.5, stop = 2.5, length = 2)
    ys = range(-2.5, stop = 2.5, length = 2)

    
    #contour!( foldl(+, interations), levels=-1:0.1:1, labels=true )
   
    function es(x, y)
        return reduce( (state, t) -> Point2f( state[1] + t[1], state[2] + t[2]), 
                                    map( q -> gradient_p(x, y, q[1], q[2], q[3]), qs), init = Point2f(0, 0))
    end   
   #=  println("Aqui estoy")
    zs = [es(x, y) for x in xs, y in ys] =#
    #println(zs)
    #arrows!( ax1, xs, ys, x -> es(x[1], x[2]), levels=-2:0.1:2, arrowsize = 7, lengthscale = 0.25 )     
    figure, ax, pl = streamplot(es, -1.2 .. 2.5, -2.5 .. 2.5, 
                                colormap = :plasma, 
                                gridsize= (56,56), 
                                arrow_size = 10,
                     )     

    #fels, sels = zip(interactions...)
    # https://beautiful.makie.org/dev/examples/generated/2d/contour/qubit/
    #fig, ax, pl = streamplot(f, -1.5..1.5, -1.5..1.5, colormap = :magma)
    # you can also pass a function to `color`, to either return a number or color value
    #streamplot(fig[1,2], f, -1.5 .. 1.5, -1.5 .. 1.5, color=(p)-> RGBAf(p..., 0.0, 1))
    #fig

   #xs = range(-15, stop = 15, length = 1000)
    #ys = range(-15, stop = 15, length = 1000)
    #U = map(q -> dp_dx(xs, ys, q[1], q[2],[3]), qs)
    #V = map(q -> dp_dy(xs, ys, q[1], q[2], q[3]), qs)
    #fig, ax, pl = streamplot(ff, -1.5..1.5, -1.5..1.5, colormap = :magma)
    #streamplot!(ax2, ff, -1.5 .. 1.5, -1.5 .. 1.5, color = :blue)
    #for q in qs
        #streamplot!(fig[1,1], (x, y) -> gradient_p(x, y, q[1], q[2], q[3]), -1.5 .. 1.5, -1.5 .. 1.5, colormap = :magma)
    #end
    uuid = UUIDs.uuid4()
    save("$uuid.png", figure)
    display(figure)
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

# Has isue, ploting in reverse
function plotting_non_numerial_aproximation()
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


#= foreach((quadrant, index) -> begin
        x_values = -9.0:00.1:9.0
        y_values = map( x -> f(x, a_values[index], b_values[index], h_values[index]), x_values)
        scatter!(
            b_values[index], 
            f(b_values[index], a_values[index], b_values[index], h_values[index]), 
            color = :red, 
            markersize = 20, 
            label = "(b, 0)"
        )
        lines!(quadrant[index], x_values, y_values, label = "cosh")
        text!(b_values[index] + 0.5, 20, text = L"(%$a,%$h)", color = :black)
        text!(-4.5, 500, text =  L"f(x)=cosh(\frac{x-b}{a}) + h", color = :black)
        vlines!(quadrant[index], [b_values[index]], color = :red, linestyle = :dash)
        hlines!(quadrant[index], [h_values[index]], xmax = [1], color = :blue, linestyle = :dash)
        display(fig)
        uuid = UUIDs.uuid4()
        save("$uuid.png", fig)    
    end, quadrants) =#


    function plot_pendulum_system()
        # Define the parameters
        g = 9.8  # Acceleration due to gravity (m/s^2)
        l = 0.5  # Length of the pendulum (m)
        gamma = 0.5  # Damping coefficient (kg/s)
        m = 1.0  # Mass of the pendulum (kg)

        # Define the function for the system of differential equations
        function pendulum(t, y)
            theta = y[1]
            theta_dot = y[2]
            return [theta_dot, -g/l * theta - gamma/m * theta_dot]
        end

        # Define the initial conditions
        theta0 = pi/5  # Initial angle (rad)
        theta_dot0 = 0.0  # Initial angular velocity (rad/s)
        y0 = [theta0, theta_dot0]

        # Define the time span
        t_start = 0.0
        t_end = 10.0
        dt = 0.01
        t = t_start:dt:t_end

        # Use the Runge-Kutta method to solve the system of differential equations
        function runge_kutta(f, t, y0)
            y = zeros(length(t), length(y0))
            y[1, :] = y0
            for i in 1:length(t)-1
                h = t[i+1] - t[i]
                k1 = f(t[i], y[i, :])
                k2 = f(t[i] + h/2, y[i, :] + h/2 * k1)
                k3 = f(t[i] + h/2, y[i, :] + h/2 * k2)
                k4 = f(t[i] + h, y[i, :] + h * k3)
                y[i+1, :] = y[i, :] + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
            end
            return y
        end

        # Solve the system of differential equations
        y = runge_kutta(pendulum, t, y0)

        # Extract the angles and angular velocities
        theta = y[:, 1]
        theta_dot = y[:, 2]
        fig = Figure()
        ax = Axis(fig[1, 1])
        lines!(ax, t, theta, label="y(x)")
        uuid = UUIDs.uuid4()
        save("$uuid.png", fig)
        display(fig)
    end    


    function plot_pendulum_system_differential_equations_library()

        # Define the parameters
        g = 9.8  # Acceleration due to gravity (m/s^2)
        l = 1.0  # Length of the pendulum (m)
        gamma = 0.5  # Damping coefficient (kg/s)
        m = 1.0  # Mass of the pendulum (kg)

        # Define the initial conditions
        theta0 = pi/5  # Initial angle (rad)
        theta_dot0 = 0.0  # Initial angular velocity (rad/s)

         # Define the function for the system of differential equations
        pendulum!( du, u, p, t) = -(g/l)sin(u) - (gamma/m)du
        # Define the time span
        t_start = 0.0
        t_end = 10.0
        tspan = (t_start, t_end)

        # Solve the system of differential equations
        prob = SecondOrderODEProblem(pendulum!, theta_dot0, theta0, tspan)
        sol = solve(prob, DPRKN6(), saveat=0.01) # Runge-Kutta-Nyström 6 order numerical method.

        # Extract the angles and angular velocities
        theta = [odeSolutionTuple[begin] for odeSolutionTuple in sol.u]
        theta_dot = [odeSolutionTuple[end] for odeSolutionTuple in sol.u]
        println(theta)
        fig = Figure()
        ax = Axis(fig[1, 1])
       
        text!(.75, 1.25, text = L"\theta(t)", color = :black)
        lines!(ax, sol.t, theta, label=L"\theta(t)")
        lines!(ax, sol.t, theta_dot, label=L"\theta^{\prime}(t)")
        fig[1, 2] = Legend(fig, ax, "Pendulum system", framevisible = false)
        uuid = UUIDs.uuid4()
        save("$uuid.png", fig)
        display(fig) 
    end



    # Trajectories plane of the pendulum.
    # This give a x = y function wrapped in a decaing sen function, like as light wave
    # Fails and incompretion to try to plot vector lines
    function linear_light_wave()

        # Define the parameters
        g = 9.8  # Acceleration due to gravity (m/s^2)
        l = 1.0  # Length of the pendulum (m)
        gamma = 0.5  # Damping coefficient (kg/s)
        m = 1.0  # Mass of the pendulum (kg)
        theta0 = pi/5  # Initial angle (rad)
        theta_dot0 = 0.0 
        # Define the vector field function for the system of differential equations
        function vectorfield!(du, u, p, t)
            theta = u[1]
            theta_dot = u[2]
            du[1] = theta_dot
            du[2] = -g/l * sin(theta) - gamma/m * theta_dot
        end 
        
        # Define the range of theta and theta_dot values
        theta_range = 0:0.1:10
        theta_dot_range = 0:0.1:10


        pendulum!( du, u, p, t) = -(g/l)sin(u) - (gamma/m)du
        # Define the time span
        t_start = -5.0
        t_end = 5.0
        tspan = (t_start, t_end)

        # Solve the system of differential equations
        prob = SecondOrderODEProblem(pendulum!, theta_dot0, theta0, tspan)
        sol = solve(prob, DPRKN6(), saveat=0.05) # Runge-Kutta-Nyström 6 order numerical method.

        # Extract the angles and angular velocities
        theta = [odeSolutionTuple[begin] for odeSolutionTuple in sol.u]
        theta_dot = [odeSolutionTuple[end] for odeSolutionTuple in sol.u]

        function meshgrid(x, y)
            X = [x for _ in y, x in x]
            Y = [y for y in y, _ in x]
            X, Y
         end
        
        # Create a meshgrid of theta and theta_dot values
        theta_mesh, theta_dot_mesh = meshgrid(theta_range, theta_dot_range)
        
        # Initialize the vector field values
        du_mesh = similar(theta_mesh)
        println(theta_mesh)

        x = 1:3
        y = 1:5
        x' .* ones(5)

        #=
        5×3 Array{Float64,2}:
        1.0  2.0  3.0
        1.0  2.0  3.0
        1.0  2.0  3.0
        1.0  2.0  3.0
        1.0  2.0  3.0
        =#
        
        ones(3)' .* y

        #=
        5×3 Array{Float64,2}:
        1.0  1.0  1.0
        2.0  2.0  2.0
        3.0  3.0  3.0
        4.0  4.0  4.0
        5.0  5.0  5.0
        =#
       #=  for i in eachindex(theta_mesh, 1), j in eachindex(theta_mesh, 2)
            vectorfield!(du_mesh[i, j, :], [theta_mesh[i, j], theta_dot_mesh[i, j]], [], 0)
        end =#
        
        # Create the scene
        println()
        # Add the quiver plot of the vector field
        us = [0.5 for x  in sol.t, _ in sol.t]
        vs = [-(g/l)sin(u) - (gamma/m)du for u in theta, du in theta_dot]
        println(vs)
        f = Figure(size = (800, 800))
        Axis(f[1, 1], backgroundcolor = "black")
        strength = vec(sqrt.(theta .^ 2 .+ theta_dot .^ 2))

        arrows!(sol.t, sol.t, theta, theta_dot, arrowsize = 10, lengthscale = 0.3,
            arrowcolor = strength, linecolor = strength)
        
        # Add the scatter plot for the critical point
   #=      scatter!(
            scene,
            [0],
            [0],
            zeros(1),
            markersize = 5,
            color = :red,
            marker = :circle,
            label = "Critical Point (theta=0, theta=0)",
        ) =#
        
        # Set the axes labels and title
  
  
        
        # Show the scene
        display(f)
    end

    function exampe_direction_field_with_Sen_cos_cuadrants_from_makie_example()
        f = Figure(size = (800, 800))
        Axis(f[1, 1], backgroundcolor = "black")

        xs = LinRange(0, 2pi, 20)
        ys = LinRange(0, 3pi, 20)
        us = [sin(x) * cos(y) for x in xs, y in ys]
        vs = [-cos(x) * sin(y) for x in xs, y in ys]
        strength = vec(sqrt.(us .^ 2 .+ vs .^ 2))

        arrows!(xs, ys, us, vs, arrowsize = 10, lengthscale = 0.3,
            arrowcolor = strength, linecolor = strength)

        display(f)
        
    end

    function ugly_pick_to_try_resolve_the_problem()
    
        # Define the differential equation
        function diff_eq(du, u, p, x)
            gamma, m, g, l, c = p
            du[1] = (u[1] / (-gamma/m * u[1] - g/l * sin(x))) 
        end
        
        # Define the parameter values
        gamma = 0.5  # damping coefficient
        m = 1.0  # mass
        g = 9.81  # acceleration due to gravity
        l = 1.0  # length of the pendulum
        
        # Define the initial condition and range of x values
        u0 = [.5]  # initial condition for y
        tspan = (-20.0,20.0)
        
        f = Figure(size = (800, 800))
        Axis(f[1, 1])

        for c in -5:1:5
             # Solve the differential equation
            prob = ODEProblem(diff_eq, u0, tspan, [gamma, m, g, l, c])
            sol = solve(prob, Tsit5(), saveat = 0.05)
            
            # Plot the solution using Makie  
            vs = collect(Iterators.flatten(sol.u))
            lines!(sol.t, vs.+c, color = :blue, linewidth = 2, label = "y(x)")
        end
       
        
        display(f)
    end    



    function non_linear_pendulum_system_direction_field()

        # Define the parameter values
        gamma = 0.5  # damping coefficient
        m = 1.0  # mass
        g = 9.81  # acceleration due to gravity
        l = 1.0  # length of the pendulum

        nonStablePoincare(x,y) = Point2f(x*(x^2+y^2-1) - y*(x^2+y^2+1),y*(x^2+y^2-1) + x*(x^2+y^2+1))
        stableVanDerPaul(x,y) = Point2f(y, -gamma/m * y - g/l * sin(x))
        semiStable(x,y) = Point2f(-y+ x*(-1+x^2+y^2)^2, x+y*(-1+x^2+y^2)^2) 

        scene = Figure(size = (800, 800))
        ax1 = Axis(scene[1, 1])
        #ax2 = Axis(scene[1, 2], backgroundcolor = "black")
        #ax3 = Axis(scene[1, 3], backgroundcolor = "black")
        
        streamplot!(stableVanDerPaul, -15..15, -15..15, colormap = :plasma, 
            gridsize= (128,128), arrow_size = 6)
        lines!( -15..15, x -> (1.5)g/l * sin(x + pi), label=L"\theta(t)", color=:orangered, linestyle = :dashdot, linewidth = 2)
        #streamplot!(ax2, stableVanDerPaul, -4..4, -4..4, colormap = :viridis, 
            #gridsize= (32,32), arrow_size = 0.07)
        #streamplot!(ax3, semiStable, -4..4, -4..4, colormap = :rainbow1, 
            #gridsize= (32,32), arrow_size = 0.07)
        #hideydecorations!(ax2, grid = false)
        #hideydecorations!(ax3, grid = false)
        limits!(ax1, -15, 15, -15, 15)
        #limits!(ax2, -4, 4, -4, 4)
        #limits!(ax3, -4, 4, -4, 4)
        display(scene)
        save("odeField3.png", scene)
        
    end

    function plot_linear_pendulum_system_with_diff()

        # Define the parameters
        g = 9.8  # Acceleration due to gravity (m/s^2)
        l = 1.0  # Length of the pendulum (m)
        gamma = 0.5  # Damping coefficient (kg/s)
        m = 1.0  # Mass of the pendulum (kg)

        # Define the initial conditions
        theta0 = pi/2  # Initial angle (rad)
        theta_dot0 = 0.0  # Initial angular velocity (rad/s)

         # Define the function for the system of differential equations
        nonlinear_pendulum!( du, u, p, t) = -(g/l)sin(u) - (gamma/m)du
        linear_pendulum!( du, u, p, t) = -(g/l)u - (gamma/m)du
        # Define the time span
        t_start = 0.0
        t_end = 10.0
        tspan = (t_start, t_end)

        # Solve the system of differential equations
        prob = SecondOrderODEProblem(nonlinear_pendulum!, theta_dot0, theta0, tspan)
        sol = solve(prob, DPRKN6(), saveat=0.01) # Runge-Kutta-Nyström 6 order numerical method.

        linear_prob = SecondOrderODEProblem(linear_pendulum!, theta_dot0, theta0, tspan)
        linear_sol = solve(linear_prob, DPRKN6(), saveat=0.01) # Runge-Kutta-Nyström 6 order numerical method.


        # Extract the angles and angular velocities
        theta = [odeSolutionTuple[begin] for odeSolutionTuple in sol.u]
        theta_dot = [odeSolutionTuple[end] for odeSolutionTuple in sol.u]

        linear_theta = [odeSolutionTuple[begin] for odeSolutionTuple in linear_sol.u]
        linear_theta_dot = [odeSolutionTuple[end] for odeSolutionTuple in linear_sol.u]


        println(theta)
        fig = Figure()
        ax = Axis(fig[1, 1])

        function diff_linear_vs_non_linear()  
            theta .- linear_theta
        end         
        text!(.75, 1.25, text = L"\theta(t)", color = :black)
        lines!(ax, sol.t, theta, label=L"non~linear~\theta(t)")
        #lines!(ax, sol.t, theta_dot, label=L"\theta^{\prime}(t)")
        lines!(ax, linear_sol.t, linear_theta, label=L"linear~\theta^{\prime}(t)")
        #lines!(ax, 0:0.01:1, x -> (g/l)*x, label=L"non~linear~f(t)=t")
        lines!(ax, linear_sol.t, diff_linear_non_linear(), label=L"diff")
        fig[1, 2] = Legend(fig, ax, "Pendulum system", framevisible = false)
        uuid = UUIDs.uuid4()
        save("$uuid.png", fig)
        display(fig) 
    end

    function plot_linear_pendulum_system_absolute_diff()

        # Define the parameters
        g = 9.8  # Acceleration due to gravity (m/s^2)
        l = 1.0  # Length of the pendulum (m)
        gamma = 0.5  # Damping coefficient (kg/s)
        m = 1.0  # Mass of the pendulum (kg)

        # Define the initial conditions
        theta0 = pi/20  # Initial angle (rad)
        theta_dot0 = 0.0  # Initial angular velocity (rad/s)

         # Define the function for the system of differential equations
        nonlinear_pendulum!( du, u, p, t) = -(g/l)sin(u) - (gamma/m)du
        linear_pendulum!( du, u, p, t) = -(g/l)u - (gamma/m)du
        # Define the time span
        t_start = 0.0
        t_end = 10.0
        tspan = (t_start, t_end)

        # Solve the system of differential equations
        prob = SecondOrderODEProblem(nonlinear_pendulum!, theta_dot0, theta0, tspan)
        sol = solve(prob, DPRKN6(), saveat=0.01) # Runge-Kutta-Nyström 6 order numerical method.

        linear_prob = SecondOrderODEProblem(linear_pendulum!, theta_dot0, theta0, tspan)
        linear_sol = solve(linear_prob, DPRKN6(), saveat=0.01) # Runge-Kutta-Nyström 6 order numerical method.


        # Extract the angles and angular velocities
        theta = [odeSolutionTuple[begin] for odeSolutionTuple in sol.u]
        theta_dot = [odeSolutionTuple[end] for odeSolutionTuple in sol.u]

        linear_theta = [odeSolutionTuple[begin] for odeSolutionTuple in linear_sol.u]
        linear_theta_dot = [odeSolutionTuple[end] for odeSolutionTuple in linear_sol.u]


        println(theta)
        fig = Figure()
        ax = Axis(fig[1, 1])

        function absolute_diff()  
            abs.(theta .- linear_theta)
        end         
        text!(.75, 1.25, text = L"\theta(t)", color = :black)
        lines!(ax, sol.t, theta, label=L"non~linear~\theta(t)")
        #lines!(ax, sol.t, theta_dot, label=L"\theta^{\prime}(t)")
        lines!(ax, linear_sol.t, linear_theta, label=L"linear~\theta^{\prime}(t)")
        #lines!(ax, 0:0.01:1, x -> (g/l)*x, label=L"non~linear~f(t)=t")
        lines!(ax, linear_sol.t, absolute_diff(), label=L"diff")
        fig[1, 2] = Legend(fig, ax, "Pendulum system", framevisible = false)
        uuid = UUIDs.uuid4()
        save("$uuid.png", fig)
        display(fig) 
    end

end # module edo_task_2   


#= obra "Física", en la cual abordó el tema del movimiento desde una perspectiva filosófica y conceptual. Su enfoque en la mecánica se diferencia significativamente de las ideas modernas, ya que Aristóteles no utilizó el enfoque matemático que caracteriza a la física clásica y moderna.

Clasificación de Movimientos:
Aristóteles clasificó el movimiento en varias categorías, incluyendo:

Movimiento Natural: El movimiento intrínseco a los objetos debido a su naturaleza interna. Por ejemplo, una piedra cae hacia abajo debido a su naturaleza pesada.
Movimiento Violento: El movimiento impuesto a un objeto por una fuerza externa. Por ejemplo, empujar una piedra.

Tipos de Movimiento:
Aristóteles también clasificó el movimiento en dos tipos:

Movimiento Natural: Aquel que ocurre de acuerdo con la naturaleza inherente de los objetos. Por ejemplo, un objeto pesado cae hacia abajo debido a su naturaleza gravitatoria.
Movimiento Violento: Aquel que se produce por una fuerza externa y va en contra de la naturaleza intrínseca de un objeto.
Lugar y Movimiento:
Aristóteles también abordó la relación entre el lugar y el movimiento. Creía que cada objeto tiene un lugar natural al que tiende a moverse.

Es importante destacar que la filosofía de Aristóteles sobre el movimiento y el cambio fue influencial en su tiempo, pero con el advenimiento de la física moderna, especialmente a partir de las ideas de Galileo y Newton, muchas de sus concepciones fueron revisadas y modificadas. La física clásica proporcionó una nueva comprensión matemática del movimiento y el cambio que superó las ideas aristotélicas en muchos aspectos.




User
mecanica del movimiento en aristoteles
ChatGPT
Aristóteles desarrolló una teoría de la mecánica del movimiento en su obra "Física", en la cual abordó el tema del movimiento desde una perspectiva filosófica y conceptual. Su enfoque en la mecánica se diferencia significativamente de las ideas modernas, ya que Aristóteles no utilizó el enfoque matemático que caracteriza a la física clásica y moderna.

A continuación, se presentan algunos puntos clave de la mecánica del movimiento según Aristóteles:

Clasificación de Movimientos:
Aristóteles clasificó el movimiento en varias categorías, incluyendo:

Movimiento Natural: El movimiento intrínseco a los objetos debido a su naturaleza interna. Por ejemplo, una piedra cae hacia abajo debido a su naturaleza pesada.
Movimiento Violento: El movimiento impuesto a un objeto por una fuerza externa. Por ejemplo, empujar una piedra.
Lugar Natural y Lugar Forzado:
Aristóteles asoció cada objeto con un "lugar natural" hacia el cual tiende a moverse. Por ejemplo, el lugar natural de un objeto pesado es hacia abajo. Además, un objeto podría ser movido a un "lugar forzado" por una fuerza externa.

Causas del Movimiento:
Aristóteles identificó cuatro causas para explicar el movimiento:

Causa Material: La sustancia o materia de la cual está hecho el objeto.
Causa Formal: La forma o estructura del objeto.
Causa Eficiente: La fuente que causa el movimiento.
Causa Final: El propósito o el fin del movimiento.
Velocidad y Aceleración:
A diferencia de las teorías modernas, Aristóteles no desarrolló una noción cuantitativa de velocidad y aceleración. Su enfoque estaba más orientado hacia las cualidades generales del movimiento que hacia mediciones precisas.

Aristóteles tenía una teoría del movimiento natural que formaba parte de su filosofía física. Su enfoque se basaba en la observación y clasificación de los cambios en el mundo natural. Aquí hay algunos puntos clave sobre el movimiento natural en la filosofía aristotélica:

Tipos de Movimiento:

Movimiento Natural: Era aquel que ocurría de forma inherente a los objetos debido a su naturaleza interna. Aristóteles creía que cada tipo de materia tenía una tendencia intrínseca a moverse de una manera específica. Por ejemplo, un objeto pesado caería hacia abajo porque su naturaleza era buscar su lugar natural en la Tierra.
Movimiento Violento: Era el movimiento que se producía debido a una fuerza externa que actuaba sobre un objeto y que iba en contra de su naturaleza interna. Por ejemplo, lanzar una piedra hacia arriba.
Lugar Natural:

Aristóteles asociaba cada tipo de materia con un "lugar natural" hacia el cual tendía a moverse. Por ejemplo, el lugar natural del fuego sería hacia arriba, mientras que el lugar natural de la tierra sería hacia abajo.
Causas del Movimiento:

Aristóteles identificó cuatro causas para explicar el movimiento:
Causa Material: La sustancia o materia de la cual está hecho el objeto.
Causa Formal: La forma o estructura del objeto.
Causa Eficiente: La fuente que causa el movimiento.
Causa Final: El propósito o el fin del movimiento.
Teleología:

La teoría aristotélica del movimiento estaba imbuida de teleología, que es la idea de que los objetos naturales tienen un propósito o finalidad intrínseca. Por ejemplo, un objeto pesado cae hacia abajo para alcanzar su lugar natural.
Es importante señalar que las ideas de Aristóteles sobre el movimiento natural fueron aceptadas y dominantes durante mucho tiempo en la historia, pero con el advenimiento de la física moderna, especialmente las contribuciones de Galileo y Newton, estas ideas fueron revisadas y reemplazadas por un enfoque más cuantitativo y matemático del movimiento.

La física moderna es el conjunto de teorías y conceptos desarrollados desde el siglo XVII en adelante

La física moderna es el conjunto de teorías y conceptos desarrollados desde el siglo XVII en adelante, que marcaron un cambio significativo en la comprensión del universo. Uno de los primeros y más influyentes científicos en este periodo fue Galileo Galilei. A continuación, se destacan algunas de las contribuciones más importantes de Galileo a la física moderna:

Ley del Movimiento Uniforme:

Galileo desarrolló la ley del movimiento uniforme, que establece que un objeto en movimiento constante seguirá moviéndose a una velocidad constante en línea recta, a menos que una fuerza externa actúe sobre él. Esta idea sentó las bases para el principio de inercia, posteriormente desarrollado por Newton.
Experimentos sobre la Caída de los Cuerpos:

A través de experimentos cuidadosos, Galileo demostró que los objetos caen al mismo ritmo independientemente de su masa. Este descubrimiento contradecía la visión aristotélica de que los objetos más pesados caen más rápido que los objetos más ligeros.
Ley de la Inercia:

Galileo formuló la ley de la inercia, que establece que un objeto en reposo permanecerá en reposo y un objeto en movimiento continuará moviéndose a una velocidad constante en línea recta, a menos que una fuerza externa actúe sobre él. Esta ley es una de las tres leyes del movimiento de Newton.

René Descartes, filósofo, matemático y científico del siglo XVII, hizo importantes contribuciones a la mecánica y al estudio del movimiento. A continuación se presentan algunas de las ideas clave de Descartes en este ámbito:

Mecanicismo y Naturaleza como Máquina:

Descartes adoptó un enfoque mecanicista para comprender la naturaleza. Veía el universo como una máquina gigante, compuesta por partículas materiales en movimiento.
Teoría del Movimiento:

Descartes propuso una teoría del movimiento basada en la colisión de partículas materiales. Según su visión, el movimiento de los cuerpos se debe a la transferencia de movimiento de una partícula a otra mediante colisiones.
Conservación de la Cantidad de Movimiento:

Descartes formuló la ley de conservación de la cantidad de movimiento. Según esta ley, la cantidad total de movimiento en un sistema aislado se mantiene constante antes y después de las colisiones entre partículas.

Isaac Newton y Gottfried Wilhelm Leibniz fueron dos de los científicos más influyentes en el desarrollo de la física y las matemáticas durante el siglo XVII. Ambos hicieron contribuciones fundamentales al estudio del movimiento y la formulación del cálculo. A continuación, se presentan algunas de sus ideas clave relacionadas con el movimiento:

Isaac Newton:
Leyes del Movimiento:

Newton formuló las tres leyes del movimiento que son fundamentales en la física clásica:
Primera Ley (Ley de la Inercia): Un objeto permanece en reposo o en movimiento rectilíneo uniforme a menos que una fuerza externa actúe sobre él.
Segunda Ley (Fuerza y Aceleración): La fuerza neta aplicada a un objeto es igual a la masa del objeto multiplicada por su aceleración (F = m * a).
Tercera Ley (Acción y Reacción): Por cada acción, hay una reacción igual y opuesta.
Ley de la Gravitación Universal:

Newton formuló la ley de la gravitación universal, que describe la atracción entre dos objetos debido a sus masas y la distancia entre ellos.
Desarrollo del Cálculo:

Newton contribuyó al desarrollo del cálculo, especialmente en la formulación de métodos para trabajar con derivadas e integrales. Su trabajo se centró en problemas de movimiento y cambio.
Gottfried Wilhelm Leibniz:
Desarrollo Independiente del Cálculo:

Leibniz también contribuyó al desarrollo del cálculo de manera independiente de Newton. Introdujo notaciones y conceptos fundamentales, como el uso de la notación "dy/dx" para representar derivadas.
Principio de la Conservación de la Energía:

Leibniz propuso la idea del principio de la conservación de la energía, que sostiene que la cantidad total de energía en un sistema aislado se mantiene constante.
Teoría Dinámica de la Materia:

Leibniz desarrolló una teoría dinámica de la materia en la que consideraba que las partículas fundamentales poseían una especie de "fuerza viva", que luego influyó en el concepto de energía cinética.
Concepto de Mónadas:

Leibniz introdujo el concepto de mónadas, entidades indivisibles que constituían la realidad fundamental. Este concepto influyó en su visión de la dinámica y la naturaleza. =#