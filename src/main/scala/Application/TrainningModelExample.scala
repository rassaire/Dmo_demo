package Application

import java.io.File

import DmoGpm.JointTriangleMesh3D
import DmoGpm.dataset.DataCollectionOfMultiMesh
import scalismo.common.UnstructuredPointsDomain
import scalismo.geometry.{EuclideanVector3D, _3D}
import scalismo.io.{LandmarkIO, MeshIO}
import scalismo.registration.TranslationTransform
import scalismo.ui.api.ScalismoUI
import DmoGpm.io.StatisticalModelIO
import breeze.linalg.DenseVector
import scalismo.utils.Random.implicits._


object TrainningModelExample {
  def main(args: Array[String]): Unit = {
    scalismo.initialize()

    /** This example is for the lollipop joint generated, where the second and the third object are the same
      * */




    /** provide the path of the folders containing the first, the second and the third objects
      *  constituting the joint respectively  */

        val st1: String = "data\\lollipop\\data shape-kinematic with no correlation\\First objects\\"
        val st2: String = "data\\lollipop\\data shape-kinematic with no correlation\\second objects\\"
        val st3: String = "data\\lollipop\\data shape-kinematic with no correlation\\second objects\\"

       /** provide the path of the folder containing the rotation centre 1, 2 and 3 (.json file) repectively for the
         * first, the second and the third objects*/
        val rotcenter1= LandmarkIO.readLandmarksJson[_3D](new File("data\\lollipop\\rotation centres\\1_30_0.json")).get
        val rotcenter2= LandmarkIO.readLandmarksJson[_3D](new File("data\\lollipop\\rotation centres\\1_30_0.json")).get
        val rotcenter3= LandmarkIO.readLandmarksJson[_3D](new File("data\\lollipop\\rotation centres\\1_30_0.json")).get

    /** provide the path of the folder  containing the points of the reference object  that represent  pose features 1, 2
      * and 3 repectively for the
      * first, the second and the third objects: It is advise to have these point
      * from the middle (in term of length) of the object up to the joint contact
      * part of the object, a lot points is not needed*/

        val posefeature1 = LandmarkIO.readLandmarksJson[_3D](new File("data\\lollipop\\poserefencelandmarks.json")).get
        val posefeature2 = LandmarkIO.readLandmarksJson[_3D](new File("data\\lollipop\\poserefencelandmarks.json")).get
        val posefeature3 = LandmarkIO.readLandmarksJson[_3D](new File("data\\lollipop\\poserefencelandmarks.json")).get




    /** laod the reference 1, 2 and 3 repctively the refrence of the object 1, 2, and 3*/
        val referenceMesh1 = MeshIO.readMesh(new File("data\\lollipop\\data shape-kinematic with no correlation\\First objects\\Synth1_30_0.stl")).get
        val referenceMesh2 = MeshIO.readMesh(new File("data\\lollipop\\data shape-kinematic with no correlation\\second objects\\Synth1_1_0.stl")).get
        val referenceMesh3 = MeshIO.readMesh(new File("data\\lollipop\\data shape-kinematic with no correlation\\second objects\\Synth1_1_0.stl")).get







        // create a joint of three structure: actually the joint must have exactly 3 structures


        val referencejoint = JointTriangleMesh3D(3,UnstructuredPointsDomain[_3D](referenceMesh1.pointSet.points.toIndexedSeq ++ referenceMesh2.pointSet.points.toIndexedSeq ++ referenceMesh3.pointSet.points.toIndexedSeq),
          List(referenceMesh1, referenceMesh2, referenceMesh3),
          List(posefeature1.map(l=>referenceMesh1.pointSet.findClosestPoint(l.point).id).toIndexedSeq,posefeature2.map(l=>referenceMesh2.pointSet.findClosestPoint(l.point).id).toIndexedSeq
            ,posefeature3.map(l=>referenceMesh3.pointSet.findClosestPoint(l.point).id).toIndexedSeq),
          List(rotcenter1.head.point,rotcenter2.head.point,rotcenter3.head.point))


        //create a sequence of registered muli meshes
        val regitseredjointmeshes = for (i <- 0 to new File(st1).listFiles().size - 1) yield {
          val mesh1 = MeshIO.readMesh(new File(st1).listFiles().sortBy(_.getName).apply(i)).get
          val mesh2 =MeshIO.readMesh(new File(st2).listFiles().sortBy(_.getName).apply(i)).get
          val mesh3 = MeshIO.readMesh(new File(st3).listFiles().sortBy(_.getName).apply(i)).get

          JointTriangleMesh3D(3,
          UnstructuredPointsDomain[_3D](mesh1.pointSet.points.toIndexedSeq ++ mesh2.pointSet.points.toIndexedSeq ++ mesh3.pointSet.points.toIndexedSeq),
            List(mesh1, mesh2, mesh3),List(posefeature1.map(l=>referenceMesh1.pointSet.findClosestPoint(l.point).id).toIndexedSeq,posefeature2.map(l=>referenceMesh2.pointSet.findClosestPoint(l.point).id).toIndexedSeq
              ,posefeature3.map(l=>referenceMesh3.pointSet.findClosestPoint(l.point).id).toIndexedSeq),
            List(rotcenter1.head.point,rotcenter2.head.point,rotcenter3.head.point))
        }

        //create multi object  mesh datacollection
        val dc = DataCollectionOfMultiMesh.fromMeshSequence(referencejoint, regitseredjointmeshes)._1.get

        //create model
        val model = DmoGpm.Model.createUsingPCA(dc).get

        // save the model
        StatisticalModelIO.writeModel(model, new File("data\\lollipop\\JointModel.h5"))

  }
}