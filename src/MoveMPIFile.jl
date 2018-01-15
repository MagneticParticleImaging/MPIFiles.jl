export MoveMPIFile

#=
 TODO: - copy most of MultiMPIFile
       - The type should get
         1. the numberOfPeriodsPerPatch
         2. the periods at which these patches should be taken (initial periods
            of each patch)
         3. the new acqOffsetField (one has to take the gradient into account)

=#


type MoveMPIFile <: MPIFile
  file::MPIFile
  # TODO
  # TODO...

  function MoveMPIFile(filename::String, )
    return new([MPIFile(f) for f in filenames])
  end

end
