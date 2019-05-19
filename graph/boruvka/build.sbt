scalaVersion := "2.11.12"

sparkVersion := "2.4.0"
sparkComponents ++= Seq("sql")
sparkComponents ++= Seq("graphx")

// set the main class for packaging the main jar
mainClass in (Compile, packageBin) := Some("com.wrapper.BoruvkaAlgorithm")

// set the main class for the main 'sbt run' task
mainClass in (Compile, run) := Some("com.wrapper.BoruvkaAlgorithm")