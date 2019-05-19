package com.wrapper

import org.apache.spark.sql.expressions.MutableAggregationBuffer
import org.apache.spark.sql.expressions.UserDefinedAggregateFunction
import org.apache.spark.sql.Row
import org.apache.spark.sql.types._

class CustomMin extends UserDefinedAggregateFunction {
  // This is the input fields for your aggregate function.
  override def inputSchema: org.apache.spark.sql.types.StructType = StructType(
        StructField("dst", IntegerType, true) ::
        StructField("weight", DoubleType, true) ::
        StructField("UniqueID", LongType, true) ::  Nil
  )

  // This is the internal fields you keep for computing your aggregate.
  override def bufferSchema: StructType = StructType(
        StructField("dstMin", IntegerType, true) ::
        StructField("weightMin", DoubleType, true) ::
        StructField("UniqueIDMin", LongType, true) ::  Nil
  )

  // This is the output type of your aggregatation function.
  override def dataType: DataType = ArrayType(DoubleType)
  override def deterministic: Boolean = true

  // This is the initial value for your buffer schema.
  // TODO
  override def initialize(buffer: MutableAggregationBuffer): Unit = {
    buffer(0) = 0
    buffer(1) = Double.MaxValue
    buffer(2) = 0L
  }

  // This is how to update your buffer schema given an input.
  override def update(buffer: MutableAggregationBuffer, input: Row): Unit = {
    if (buffer.getAs[Double](1) > input.getAs[Double](1)) {
        buffer(1) = input.getAs[Double](1)
        buffer(0) = input.getAs[Int](0)
        buffer(2) = input.getAs[Long](2)
    }
  }

  // This is how to merge two objects with the bufferSchema type.
  override def merge(buffer1: MutableAggregationBuffer, buffer2: Row): Unit = {
    if (buffer1.getAs[Double](1) > buffer2.getAs[Double](1)) {
        buffer1(1) = buffer2.getAs[Double](1)
        buffer1(0) = buffer2.getAs[Int](0)
        buffer1(2) = buffer2.getAs[Long](2)
    }
  }

  // This is where you output the final value, given the final value of your bufferSchema.
  override def evaluate(buffer: Row): Any = {
    Array(buffer.getAs[Int](0).toDouble, buffer.getAs[Double](1), buffer.getAs[Long](2).toDouble)
  }
}