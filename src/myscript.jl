if length(ARGS) < 2
    println("Usage: julia myscript.jl arg1 arg2")
    println("Please provide at least two command-line arguments.")
    exit(1) # Exit with error
end

# Accessing command-line arguments
arg1 = ARGS[1]
arg2 = ARGS[2]

# Print the command-line arguments
println("First argument: $arg1")
println("Second argument: $arg2")