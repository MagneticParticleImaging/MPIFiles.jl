@testset "MagneticFieldMeasurement" begin
  # Generate data
  description_ = "Just some testset."

  shp = [3, 3, 3]
  fov = [3.0, 3.0, 3.0]u"mm"
  ctr = [0.0, 0.0, 0.0]u"mm"
  positions_ = RegularGridPositions(shp, fov, ctr)

  gradient = [1.0u"T/m", 1.5u"T/m", 2.0u"T/m"]
  field_ = fill(0.0u"T", (length(positions_), 3))
  for (idx, pos) in enumerate(positions_)
    field_[idx, :] = pos.*gradient
  end
  
  fieldError_ = fill(0.1u"mT" .|> u"T", (length(positions_), 3))
  fieldFrequency_ = fill(0.0u"Hz", length(positions_))

  currents_ = fill(1.0u"A", (length(positions_), 3))
  timestamp_ = DateTime("2021-08-25T07:08:21.881")
  sensorOffset_ = [0.1u"mm", 0.2u"mm", 0.3u"mm"] .|> u"m"
  temperature_ = fill(20.0u"°C", length(positions_))

  # Set data
  measurement = MagneticFieldMeasurement()

  # Test setters
  description(measurement, description_)
  positions(measurement, positions_)
  field(measurement, field_)
  fieldError(measurement, fieldError_)
  fieldFrequency(measurement, fieldFrequency_)
  currents(measurement, currents_)
  timestamp(measurement, timestamp_)
  sensorOffset(measurement, sensorOffset_)
  temperature(measurement, temperature_)

  @test measurement.description == description_
  @test measurement.positions == positions_
  @test measurement.field == field_
  @test measurement.fieldError == fieldError_
  @test measurement.fieldFrequency == fieldFrequency_
  @test measurement.currents == currents_
  @test measurement.timestamp == timestamp_
  @test measurement.sensorOffset == sensorOffset_
  @test measurement.temperature == temperature_

  # Test getters
  @test description(measurement) == description_
  @test positions(measurement) == positions_
  @test field(measurement) == field_
  @test fieldError(measurement) == fieldError_
  @test fieldFrequency(measurement) == fieldFrequency_
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
  @test field(measurementRead) == field_
  @test fieldError(measurementRead) == fieldError_
  @test fieldFrequency(measurementRead) == fieldFrequency_
  @test currents(measurementRead) == currents_
  @test timestamp(measurementRead) == timestamp_
  @test sensorOffset(measurementRead) == sensorOffset_
  @test temperature(measurementRead) == temperature_
end