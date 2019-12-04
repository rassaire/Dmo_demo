package DmoGpm.dataset

import scalismo.common.UnstructuredPointsDomain
import scalismo.geometry.{Point, _3D}
import scalismo.mesh.TriangleMesh
import scalismo.registration.Transformation
import scalismo.mesh.TetrahedralMesh

import scala.util.{Failure, Success, Try}

private object DataUtils {
  /**
    * Partitions a list os possible transformation (tries) into those that succeeded and those who failed
    */
  def partitionSuccAndFailedTries[A](tries: Seq[Try[A]]): (Seq[A], Seq[Throwable]) = {
    val (s, f) = tries.partition(_.isSuccess)
    val throwables = for (failed <- f) yield {
      failed match {
        case Failure(t) => t
        case _ => new Throwable("This will never happen")
      }
    }
    val transforms = s.map(_.get)
    (transforms, throwables)
  }

  /**
    * Create a transformation from a mesh. The transformation maps from the reference mesh to the corresponding target point.
    */
  def meshToTransformation(refMesh: TriangleMesh[_3D], targetMesh: TriangleMesh[_3D]): Try[Transformation[_3D]] = {
    if (refMesh.pointSet.numberOfPoints != targetMesh.pointSet.numberOfPoints)
      Failure(new Throwable(s"reference and target mesh do not have the same number of points (${refMesh.pointSet.numberOfPoints} != ${targetMesh.pointSet.numberOfPoints}"))
    else {
      val t = new Transformation[_3D] {
        override val domain = refMesh.boundingBox
        override val f = (x: Point[_3D]) => {
          val ptId = refMesh.pointSet.findClosestPoint(x).id
          targetMesh.pointSet.point(ptId)
        }
      }
      Success(t)
    }
  }

  /**
    * Create a transformation from a mesh volume. The transformation maps from the reference mesh volume to the corresponding target point.
    */
  def meshVolumeToTransformation(refMesh: TetrahedralMesh[_3D], targetMesh: TetrahedralMesh[_3D]): Try[Transformation[_3D]] = {
    if (refMesh.pointSet.numberOfPoints != targetMesh.pointSet.numberOfPoints)
      Failure(new Throwable(s"reference and target mesh do not have the same number of points (${refMesh.pointSet.numberOfPoints} != ${targetMesh.pointSet.numberOfPoints}"))
    else {
      val t = new Transformation[_3D] {
        override val domain = refMesh.boundingBox
        override val f = (x: Point[_3D]) => {
          val ptId = refMesh.pointSet.findClosestPoint(x).id
          targetMesh.pointSet.point(ptId)
        }
      }
      Success(t)
    }
  }


  /**
    * Create a transformation from a mesh volume. The transformation maps from the reference mesh volume to the corresponding target point.
    */
  def multimeshToTransformation(refMesh: UnstructuredPointsDomain[_3D], targetMesh: UnstructuredPointsDomain[_3D]): Try[Transformation[_3D]] = {
    if (refMesh.numberOfPoints != targetMesh.numberOfPoints)
      Failure(new Throwable(s"reference and target mesh do not have the same number of points (${refMesh.numberOfPoints} != ${targetMesh.numberOfPoints}"))
    else {

      val t = new Transformation[_3D] {
        override val domain = refMesh.boundingBox
        override val f = (x: Point[_3D]) => {
          val ptId = refMesh.findClosestPoint(x).id
          targetMesh.point(ptId)
        }
      }
      Success(t)
    }
  }

}