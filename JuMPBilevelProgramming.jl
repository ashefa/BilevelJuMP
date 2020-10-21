#######################################################################################
# bilevel
# Solve the linear bilevel programming
#  {min/max}  dot(c1,x) + dot(d1,y)
# subject to  A1 x + B1 y <= b1
#                x >= 0
#            {min/max}  dot(c2,x) + dot(d2,y)
#           subject to  A2 x + B2 y <= b2
#                                 y >= 0
# where x and y are the vectors of upper-level and lower-level variables, respectively.
#######################################################################################

# KKT method:
function bilevel(c1,d1,c2,d2,A1,B1,A2,B2,b1,b2)
    
    m1=size(A1,1);m2=size(A2,1);nA=size(A1,2);nB=size(B1,2);
    
    using JuMP

    BiLevel = Model() 
    @variable(BiLevel, x[1:nA] ) 
    @variable(BiLevel, y[1:nB] ) 
    @variable(BiLevel, u[1:m2] ) 
    @variable(BiLevel, v[1:nB] ) 
    @variable(BiLevel, z[1:m2] , Bin ) 
    M = 1.0e8;
    @constraint(BiLevel, [A1';B1']'*[x;y] .<= b1 ) 
    @constraint(BiLevel, (u'*B2)'-v .== -d2 ) 
    @constraint(BiLevel, u .<= M*z ) 
    @constraint(BiLevel, b2-[A2';B2']'*[x;y] .<= M*(1-z) ) 
    @constraint(BiLevel, [A2';B2']'*[x;y] .<= b2 ) 
    @objective(BiLevel, Min, dot(c1',x') + dot(d1',y') ) 
    status = solve(BiLevel)

    println("x = ",getvalue(x)),println("y = ",getvalue(y)),println("obj(upper-level) = ",getobjectivevalue(BiLevel));
    
end

# Second method:
function bilevel(c1,d1,c2,d2,A1,B1,A2,B2,b1,b2)
    
    m1=size(A1,1);m2=size(A2,1);nA=size(A1,2);nB=size(B1,2);
    
    using JuMP

    BiLevel = Model() 
    @variable(BiLevel, x[1:nA] ) 
    @variable(BiLevel, y[1:nB] ) 
    @variable(BiLevel, u[1:m2] ) 
    @constraint(BiLevel, [A1';B1']'*[x;y] .<= b1 ) 
    @constraint(BiLevel, [A2';B2']'*[x;y] .<= b2 ) 
    @constraint(BiLevel, u'*[A2';B2']' .<= [c2;d2]' ) 
    @constraint(BiLevel, dot(c2',x') + dot(d2',y') == dot(b2',u') ) 
    @objective(BiLevel, Min, dot(c1',x') + dot(d1',y') ) 
    status = solve(BiLevel)
    
    println("x = ",getvalue(x)),println("y = ",getvalue(y)),println("obj(upper-level) = ",getobjectivevalue(BiLevel));
    
end