@testset "MagneticFieldMeasurement" begin
  # Generate data
  description_ = "Just some testset."

  shp = [3, 3, 3]
  fov = [3.0, 3.0, 3.0]u"mm"
  ctr = [0.0, 0.0, 0.0]u"mm"
  positions_ = RegularGridPositions(shp, fov, ctr)

  gradient = [1.0u"T/m", 1.5u"T/m", 2.0u"T/m"]
  fields_ = fill(0.0u"T", (length(positions_), 3))
  for (idx, pos) in enumerate(positions_)
    fields_[idx, :] = pos.*gradient
  end
  
  fieldsError_ = fill(0.1u"mT" .|> u"T", (length(positions_), 3))
  fieldsFrequency_ = fill(0.0u"Hz", length(positions_))

  currents_ = fill(1.0u"A", (length(positions_), 3))
  timestamp_ = DateTime("2021-08-25T07:08:21.881")
  sensorOffset_ = [0.1u"mm", 0.2u"mm", 0.3u"mm"] .|> u"m"
  temperature_ = fill(20.0u"°C", length(positions_))

  # Set data
  measurement = MagneticFieldMeasurement()

  # Test setters
  description(measurement, description_)
  positions(measurement, positions_)
  fields(measurement, fields_)
  fieldsError(measurement, fieldsError_)
  fieldsFrequency(measurement, fieldsFrequency_)
  currents(measurement, currents_)
  timestamp(measurement, timestamp_)
  sensorOffset(measurement, sensorOffset_)
  temperature(measurement, temperature_)

  @test measurement.description == description_
  @test measurement.positions == positions_
  @test measurement.fields == fields_
  @test measurement.fieldsError == fieldsError_
  @test measurement.fieldsFrequency == fieldsFrequency_
  @test measurement.currents == currents_
  @test measurement.timestamp == timestamp_
  @test measurement.sensorOffset == sensorOffset_
  @test measurement.temperature == temperature_

  # Test getters
  @test description(measurement) == description_
  @test positions(measurement) == positions_
  @test fields(measurement) == fields_
  @test fieldsError(measurement) == fieldsError_
  @test fieldsFrequency(measurement) == fieldsFrequency_
  @test currents(measurement) == currents_
  @test timestamp(measurement) == timestamp_
  @test sensorOffset(measurement) == sensorOffset_
  @test temperature(measurement) == temperature_

  # Save measurement
  magneticMeasurementBaseDir = joinpath(tmpdir, "magnetic")
  mkpath(magneticMeasurementBaseDir)
  magneticMeasurementPath = joinpath(magneticMeasurementBaseDir, "measurement.h5")
  saveMagneticFieldAsHDF5(measurement, magneticMeasurementPath)
  @test isfile(magneticMeasurementPath)

  # Read measurement back
  measurementRead = MagneticFieldMeasurement(magneticMeasurementPath)

  # Test read data
  @test description(measurementRead) == description_
  @test positions(measurementRead) == positions_
  @test fields(measurementRead) == fields_
  @test fieldsError(measurementRead) == fieldsError_
  @test fieldsFrequency(measurementRead) == fieldsFrequency_
  @test currents(measurementRead) == currents_
  @test timestamp(measurementRead) == timestamp_
  @test sensorOffset(measurementRead) == sensorOffset_
  @test temperature(measurementRead) == temperature_
end