package barneshut

import java.util.concurrent._
import scala.collection._
import scala.math._
import scala.collection.parallel._
import barneshut.conctrees.ConcBuffer
import org.junit._
import org.junit.Assert.{assertEquals, fail}

class BarnesHutSuite {
  // test cases for quad tree

import FloatOps._
  @Test def `Empty: center of mass should be the center of the cell`: Unit = {
    val quad = Empty(51f, 46.3f, 5f)
    assert(quad.massX == 51f, s"${quad.massX} should be 51f")
    assert(quad.massY == 46.3f, s"${quad.massY} should be 46.3f")
  }

  @Test def `Empty: mass should be 0`: Unit = {
    val quad = Empty(51f, 46.3f, 5f)
    assert(quad.mass == 0f, s"${quad.mass} should be 0f")
  }

  @Test def `Empty: total should be 0`: Unit = {
    val quad = Empty(51f, 46.3f, 5f)
    assert(quad.total == 0, s"${quad.total} should be 0")
  }

  @Test def `Leaf with 1 body`: Unit = {
    val b = new Body(123f, 18f, 26f, 0f, 0f)
    val quad = Leaf(17.5f, 27.5f, 5f, Seq(b))

    assert(quad.mass ~= 123f, s"${quad.mass} should be 123f")
    assert(quad.massX ~= 18f, s"${quad.massX} should be 18f")
    assert(quad.massY ~= 26f, s"${quad.massY} should be 26f")
    assert(quad.total == 1, s"${quad.total} should be 1")
  }


  @Test def `Fork with 3 empty quadrants and 1 leaf (nw)`: Unit = {
    val b = new Body(123f, 18f, 26f, 0f, 0f)
    val nw = Leaf(17.5f, 27.5f, 5f, Seq(b))
    val ne = Empty(22.5f, 27.5f, 5f)
    val sw = Empty(17.5f, 32.5f, 5f)
    val se = Empty(22.5f, 32.5f, 5f)
    val quad = Fork(nw, ne, sw, se)

    assert(quad.centerX == 20f, s"${quad.centerX} should be 20f")
    assert(quad.centerY == 30f, s"${quad.centerY} should be 30f")
    assert(quad.mass ~= 123f, s"${quad.mass} should be 123f")
    assert(quad.massX ~= 18f, s"${quad.massX} should be 18f")
    assert(quad.massY ~= 26f, s"${quad.massY} should be 26f")
    assert(quad.total == 1, s"${quad.total} should be 1")
  }

  @Test def `Empty.insert(b) should return a Leaf with only that body (2pts)`: Unit = {
    val quad = Empty(51f, 46.3f, 5f)
    val b = new Body(3f, 54f, 46f, 0f, 0f)
    val inserted = quad.insert(b)
    inserted match {
      case Leaf(centerX, centerY, size, bodies) =>
        assert(centerX == 51f, s"$centerX should be 51f")
        assert(centerY == 46.3f, s"$centerY should be 46.3f")
        assert(size == 5f, s"$size should be 5f")
        assert(bodies == Seq(b), s"$bodies should contain only the inserted body")
      case _ =>
        fail("Empty.insert() should have returned a Leaf, was $inserted")
    }
  }

  // test cases for Body

  @Test def `Body.updated should do nothing for Empty quad trees`: Unit = {
    val b1 = new Body(123f, 18f, 26f, 0f, 0f)
    val body = b1.updated(Empty(50f, 60f, 5f))

    assertEquals(0f, body.xspeed, precisionThreshold)
    assertEquals(0f, body.yspeed, precisionThreshold)
  }

  @Test def `Body.updated should take bodies in a Leaf into account (2pts)`: Unit = {
    val b1 = new Body(123f, 18f, 26f, 0f, 0f)
    val b2 = new Body(524.5f, 24.5f, 25.5f, 0f, 0f)
    val b3 = new Body(245f, 22.4f, 41f, 0f, 0f)

    val quad = Leaf(15f, 30f, 20f, Seq(b2, b3))

    val body = b1.updated(quad)

    assert(body.xspeed ~= 12.587037f)
    assert(body.yspeed ~= 0.015557117f)
  }

  // test cases for sector matrix

  @Test def `'SectorMatrix.+=' should add a body at (25,47) to the correct bucket of a sector matrix of size 96 (2pts)`: Unit = {
    val body = new Body(5, 25, 47, 0.1f, 0.1f)
    val boundaries = new Boundaries()
    boundaries.minX = 1
    boundaries.minY = 1
    boundaries.maxX = 97
    boundaries.maxY = 97
    val sm = new SectorMatrix(boundaries, SECTOR_PRECISION)
    sm += body
    val res = sm(2, 3).size == 1 && sm(2, 3).find(_ == body).isDefined
    assert(res, s"Body not found in the right sector")
  }

  @Test def `simulator updateBoundaries works`: Unit = {
    val model = new SimulationModel
    val s = new Simulator(model.taskSupport, model.timeStats)
    val boundaries = new Boundaries
    boundaries.minX = 1f
    boundaries.maxX = 5f
    boundaries.minY = 1f
    boundaries.maxY = 5f

    val bodyInside = new Body(1f, 3f, 3f, 1f, 1f)
    val bodyLeft = new Body(1f, 0f, 3f, 1f, 1f)
    val bodyRight = new Body(1f, 6f, 3f, 1f, 1f)
    val bodyAbove = new Body(1f, 3f, 0f, 1f, 1f)
    val bodyBelow = new Body(1f, 3f, 6f, 1f, 1f)

    assert(compareBoundaries(s.updateBoundaries(boundaries, bodyInside), boundaries), s"Body inside")
    val boundariesLeft = new Boundaries
    boundariesLeft.minX = 0f
    boundariesLeft.maxX = 5f
    boundariesLeft.minY = 1f
    boundariesLeft.maxY = 5f
    assert(compareBoundaries(s.updateBoundaries(boundaries, bodyLeft), boundariesLeft), s"Body left")
    val boundariesRight = new Boundaries
    boundariesRight.minX = 1f
    boundariesRight.maxX = 6f
    boundariesRight.minY = 1f
    boundariesRight.maxY = 5f
    assert(compareBoundaries(s.updateBoundaries(boundaries, bodyRight), boundariesRight), s"Body right")
    val boundariesAbove = new Boundaries
    boundariesAbove.minX = 1f
    boundariesAbove.maxX = 5f
    boundariesAbove.minY = 0f
    boundariesAbove.maxY = 5f
    assert(compareBoundaries(s.updateBoundaries(boundaries, bodyAbove), boundariesAbove), s"Body above")
    val boundariesBelow = new Boundaries
    boundariesBelow.minX = 1f
    boundariesBelow.maxX = 5f
    boundariesBelow.minY = 1f
    boundariesBelow.maxY = 6f
    assert(compareBoundaries(s.updateBoundaries(boundaries, bodyBelow), boundariesBelow), s"Body below")
  }

  @Test def `simulator mergeBoundaries works`: Unit = {
    val model = new SimulationModel
    val s = new Simulator(model.taskSupport, model.timeStats)
    val a = new Boundaries
    a.minX = 1f
    a.maxX = 6f
    a.minY = 0f
    a.maxY = 5f

    val b = new Boundaries
    b.minX = 0f
    b.maxX = 5f
    b.minY = 1f
    b.maxY = 6f

    val boundaries = new Boundaries
    boundaries.minX = 0f
    boundaries.maxX = 6f
    boundaries.minY = 0f
    boundaries.maxY = 6f
    assert(compareBoundaries(s.mergeBoundaries(a, b), boundaries), s"Body inside")
  }

  @Rule def individualTestTimeout = new org.junit.rules.Timeout(10 * 1000)

  def compareBoundaries(a: Boundaries, b: Boundaries): Boolean = {
    a.minX == b.minX && a.maxX == b.maxX && a.minY == b.minY && a.maxY == b.maxY
  }
}

object FloatOps {
  val precisionThreshold = 1e-4

  /** Floating comparison: assert(float ~= 1.7f). */
  implicit class FloatOps(val self: Float) extends AnyVal {
    def ~=(that: Float): Boolean =
      abs(self - that) < precisionThreshold
  }

  /** Long floating comparison: assert(double ~= 1.7). */
  implicit class DoubleOps(val self: Double) extends AnyVal {
    def ~=(that: Double): Boolean =
      abs(self - that) < precisionThreshold
  }

  /** Floating sequences comparison: assert(floatSeq ~= Seq(0.5f, 1.7f). */
  implicit class FloatSequenceOps(val self: Seq[Float]) extends AnyVal {
    def ~=(that: Seq[Float]): Boolean =
      self.size == that.size &&
        self.zip(that).forall { case (a, b) =>
          abs(a - b) < precisionThreshold
        }
  }
}

