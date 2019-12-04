package DmoGpm
import scalismo.common.{PointId, UnstructuredPointsDomain}
import scalismo.geometry.{NDSpace, Point, _3D}
import scalismo.mesh._
import scalismo.mesh.TetrahedralMesh
trait JointTriangleMesh[D] {
  def nomberOfobject: Int
  def pointSet: UnstructuredPointsDomain[D]
  def objects: List[TriangleMesh[D]]
  def poseIds: List[IndexedSeq[PointId]]
  def rotCenters: List[Point[_3D]]
}

trait JointTetrahedralMesh[D] {
  def nomberOfobject: Int
  def pointSet: UnstructuredPointsDomain[D]
  def objects: List[TetrahedralMesh[D]]
}

/** Standard 3D joint mesh, geometry only */
case class JointTriangleMesh3D(nomberOfobject: Int, pointSet: UnstructuredPointsDomain[_3D], objects: List[TriangleMesh[_3D]], poseIds:List[IndexedSeq[PointId]], rotCenters:List[Point[_3D]]) extends JointTriangleMesh[_3D] {

  val numberOfobject = nomberOfobject

  val domain = pointSet

  val jointmeshes = objects

  val poseids=poseIds
  val rotcenter=rotCenters

}

case class JointTetrahedralMesh3D(nomberOfobject: Int = 3, pointSet: UnstructuredPointsDomain[_3D], objects: List[TetrahedralMesh[_3D]]) extends JointTetrahedralMesh[_3D] {

  val numberOfobject = nomberOfobject

  val domain = pointSet

  val jointmeshes = objects

  lazy val boundingBox = pointSet.boundingBox
}

object JointTriangleMesh {

  def apply[D: NDSpace](objects: List[scalismo.mesh.TriangleMesh[D]],poseIds:List[IndexedSeq[PointId]], rotCenters:List[Point[_3D]])(implicit creator: Create[D]) = {
    creator.createJointTriangleMesh(objects,poseIds,rotCenters)
  }

  /** Typeclass for creating domains of arbitrary dimensionality */
  trait Create[D] extends UnstructuredPointsDomain.Create[D] {
    def createJointTriangleMesh(objects: List[scalismo.mesh.TriangleMesh[D]],poseIds:List[IndexedSeq[PointId]], rotCenters:List[Point[_3D]]): JointTriangleMesh[D]
  }

  trait Create3D extends Create[_3D] {
    override def createJointTriangleMesh(objects: List[scalismo.mesh.TriangleMesh[_3D]],poseIds:List[IndexedSeq[PointId]], rotCenters:List[Point[_3D]]) = {
      var d = List[Point[_3D]]().toIndexedSeq
      for (i <- 0 to objects.size - 1) {
        d = d ++ objects(i).pointSet.points.toIndexedSeq
      }
      JointTriangleMesh3D(objects.size, UnstructuredPointsDomain[_3D](d), objects,poseIds,rotCenters)
    }
  }

  implicit def parametricToConcreteType3D(jointMeshes: JointTriangleMesh[_3D]): JointTriangleMesh3D = {
    jointMeshes.asInstanceOf[JointTriangleMesh3D]
  }

}

object JointTetrahedralMesh {

  def apply[D: NDSpace](objects: List[TetrahedralMesh[D]])(implicit creator: Create[D]) = {
    creator.createJointTetrahedralMesh(objects)
  }

  /** Typeclass for creating domains of arbitrary dimensionality */
  trait Create[D] extends UnstructuredPointsDomain.Create[D] {
    def createJointTetrahedralMesh(objects: List[TetrahedralMesh[D]]): JointTetrahedralMesh[D]
  }

  trait Create3D extends Create[_3D] {
    override def createJointTetrahedralMesh(objects: List[TetrahedralMesh[_3D]]) = {
      var d = List[Point[_3D]]().toIndexedSeq
      for (i <- 0 to objects.size - 1) {
        d = d ++ objects(i).pointSet.points.toIndexedSeq
      }
      JointTetrahedralMesh3D(objects.size, UnstructuredPointsDomain[_3D](d), objects)
    }
  }

  implicit def parametricToConcreteType3D(jointMeshes: JointTetrahedralMesh[_3D]): JointTetrahedralMesh3D = {
    jointMeshes.asInstanceOf[JointTetrahedralMesh3D]
  }

}
