using MPIFiles
using Base.Test
using Unitful

randPos, randPosTime, numFramesToAcquire, seed = createRandPositions([10.0u"mm"; 10.0u"mm"; 10.0u"mm"], UInt64(100), UInt64(100))
randPos, randPosTime, numFramesToAcquire, seed = createRandPositions([10.0u"mm"; 10.0u"mm"; 10.0u"mm"], UInt64(100), UInt64(100),
                                                  offset =[10.0u"mm";10.0u"mm";1.0u"mm"], waitLingerTime=0.5u"s",seed=UInt64(10))
