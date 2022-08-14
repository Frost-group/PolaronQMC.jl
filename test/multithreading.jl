begin
using Base.Threads
end

begin
    nthreads()
    @threads for i in 1:10
        println([i,threadid()])
    end



end
