package DmoGpm.dataset

import java.io.File

import scalismo.common.UnstructuredPointsDomain
import scalismo.geometry._
import scalismo.io.MeshIO
import scalismo.mesh._
import scalismo.registration.{LandmarkRegistration, RigidTransformation, Transformation, TranslationTransform}
import scalismo.mesh.{TetrahedralMesh, TetrahedralMesh3D}
import scalismo.utils.Random
import DmoGpm.{JointTriangleMesh, JointTriangleMesh3D}

import scala.annotation.tailrec

case class CrossvalidationFoldMultiMesh(trainingData: DataCollectionOfMultiMesh, testingData: DataCollectionOfMultiMesh)

/**
  * A registered item in a dataset.
  *
  *  @param info A human-readable description of the processing the data item went through. Current implemented methods on data collections,
  *  such as  will increment this description
  *  @param transformation Transformation to apply to obtain the data item from the reference of the reference item of the dataset.
  *  This would typically be the transformation resulting from registering a reference mesh to the mesh represented by this data item.
  */
case class DataItem[D](info: String, transformation: Transformation[D])

/**
  * Data-structure for handling a dataset of registered 3D meshes. All pre-implemented operations such as building a
  * PCA model or performing a Generalized Procrustes Analysis require a DataCollection as input
  *
  * @param reference The reference mesh of the dataset. This is the mesh that was registered to all other items of the dataset.
  * @param dataItems Sequence of data items containing the required transformations to apply to the reference mesh in order to obtain
  * other elements of the dataset.
  */


/**
  * Data-structure for handling a dataset of registered multi 3D meshes. All pre-implemented operations such as building a
  * PCA model or performing a Generalized Procrustes Analysis require a DataCollectionofMultiMesh as input
  *
  * @param reference The reference Joint triangle  mesh of the dataset. This is the Joint triangle mesh that was registered to all other items of the dataset.
  * @param dataItems Sequence of data items containing the required transformations to apply to the Joint reference mesh in order to obtain
  * other elements of the dataset.
  */
case class DataCollectionOfMultiMesh(shapeReference: JointTriangleMesh[_3D], dataItems: Seq[DataItem[_3D]])(implicit random: Random) {
  val size: Int = dataItems.size

  val refrot1 =shapeReference.poseIds(0).map{pi=>
    val p=shapeReference.objects(0).pointSet.point(pi)
    TranslationTransform(EuclideanVector3D(-0.01, 0.0, 0.0)).apply(p)}.toIndexedSeq


  val refrot2 =shapeReference.poseIds(1).map{pi=>
    val p=shapeReference.objects(1).pointSet.point(pi)
    TranslationTransform(EuclideanVector3D(-0.01, 0.0, 0.0)).apply(p)}.toIndexedSeq

  val refrot3 =shapeReference.poseIds(2).map{pi=>
    val p=shapeReference.objects(2).pointSet.point(pi)
    TranslationTransform(EuclideanVector3D(-0.01, 0.0, 0.0)).apply(p)}.toIndexedSeq

  val reftrans1 = IndexedSeq(TranslationTransform(EuclideanVector3D(0.0, 0.0, 0.01)).apply(shapeReference.objects(0).pointSet.point(shapeReference.poseIds(0)(0))))

  val reftrans2 = IndexedSeq(TranslationTransform(EuclideanVector3D(0.0, 0.0, 0.01)).apply(shapeReference.objects(1).pointSet.point(shapeReference.poseIds(1)(0))))

  val reftrans3 = IndexedSeq(TranslationTransform(EuclideanVector3D(0.0, 0.0, 0.01)).apply(shapeReference.objects(2).pointSet.point(shapeReference.poseIds(2)(0))))

  private val d: IndexedSeq[Point[_3D]] = refrot1++ refrot2++ refrot3++ reftrans1++ reftrans2 ++ reftrans3++
    shapeReference.objects(0).pointSet.points.toIndexedSeq ++
    shapeReference.objects(1).pointSet.points.toIndexedSeq ++ shapeReference.objects(2).pointSet.points.toIndexedSeq

  val referenceShapePoseDomain = UnstructuredPointsDomain[_3D](d)

  private[dataset] def createCrossValidationFolds(nFolds: Int): Seq[CrossvalidationFoldMultiMesh] = {

    val shuffledDataItems = random.scalaRandom.shuffle(dataItems)
    val foldSize = shuffledDataItems.size / nFolds
    val dataGroups = shuffledDataItems.grouped(foldSize).toSeq

    val folds = for (currFold <- 0 until nFolds) yield {
      val testingDataItems = dataGroups(currFold)
      val testingCollection = DataCollectionOfMultiMesh(shapeReference, testingDataItems)
      val trainingDataItems = (dataGroups.slice(0, currFold).flatten ++: dataGroups.slice(currFold + 1, dataGroups.size).flatten)
      val trainingCollection = DataCollectionOfMultiMesh(shapeReference, trainingDataItems)

      CrossvalidationFoldMultiMesh(trainingCollection, testingCollection)
    }
    folds
  }

  private[dataset] def createLeaveOneOutFolds = createCrossValidationFolds(dataItems.size)

  /**
    * Returns a new DataCollectionofMeshvolume where the given function was applied to all data items
    */
  def mapItems(f: DataItem[_3D] => DataItem[_3D]): DataCollectionOfMultiMesh = {
    new DataCollectionOfMultiMesh(shapeReference, dataItems.map(f))
  }

  /**
    * Returns the mean surface computed by transforming the reference with all the transformations in the datacollection
    */
  def meanJoint: JointTriangleMesh[_3D] = {

    val currentmesh1 = shapeReference.objects.apply(0).transform(meanTransformation)
    val currentmesh2 = shapeReference.objects.apply(1).transform(meanTransformation)
    val currentmesh3 = shapeReference.objects.apply(2).transform(meanTransformation)

    val trans1p= reftrans1.map(m =>  meanTransformation.apply(m))
    val trans1=TranslationTransform(trans1p(0).toVector)

    val rot1refland = refrot1.map(m => new Landmark[_3D](shapeReference.objects.apply(0).pointSet.findClosestPoint(m).id.toString, meanTransformation.apply(m)))
    val rot1tagland = refrot1.map(m => new Landmark[_3D](shapeReference.objects.apply(0).pointSet.findClosestPoint(m).id.toString, m))
    val rot1 = LandmarkRegistration.rigid3DLandmarkRegistration(rot1tagland, rot1refland, shapeReference.rotCenters(0))

    val mesh1 = currentmesh1.transform(RigidTransformation(rot1.rotation, trans1)) //(rot1.rotation).transform(trans1.translation)



    val trans2p= reftrans2.map(m =>  meanTransformation.apply(m))
    val trans2=TranslationTransform(trans2p(0).toVector)

    val rot2refland = refrot2.map(m => new Landmark[_3D](shapeReference.objects.apply(1).pointSet.findClosestPoint(m).id.toString, meanTransformation.apply(m)))
    val rot2tagland = refrot2.map(m => new Landmark[_3D](shapeReference.objects.apply(1).pointSet.findClosestPoint(m).id.toString, m))
    val rot2 = LandmarkRegistration.rigid3DLandmarkRegistration(rot2tagland, rot2refland, shapeReference.rotCenters(1))

    val mesh2 = currentmesh1.transform(RigidTransformation(rot2.rotation, trans2)) //(rot1.rotation).transform(trans1.translation)






    val trans3p= reftrans3.map(m =>  meanTransformation.apply(m))
    val trans3=TranslationTransform(trans3p(0).toVector)

    val rot3refland = refrot3.map(m => new Landmark[_3D](shapeReference.objects.apply(2).pointSet.findClosestPoint(m).id.toString, meanTransformation.apply(m)))
    val rot3tagland = refrot3.map(m => new Landmark[_3D](shapeReference.objects.apply(2).pointSet.findClosestPoint(m).id.toString, m))
    val rot3 = LandmarkRegistration.rigid3DLandmarkRegistration(rot3tagland, rot3refland, shapeReference.rotCenters(1))

    val mesh3 = currentmesh3.transform(RigidTransformation(rot3.rotation, trans3)) //(rot1.rotation).transform(trans1.translation)

    JointTriangleMesh3D(3, UnstructuredPointsDomain(mesh1.pointSet.points.toIndexedSeq ++ mesh2.pointSet.points.toIndexedSeq ++ mesh3.pointSet.points.toIndexedSeq),
      List(mesh1, mesh2, mesh3),shapeReference.poseIds,shapeReference.rotCenters)
  }

  /**
    * Returns the mean transformation from all the transformation in the datacollectionOfMultiMesh
    */
  val meanTransformation: Transformation[_3D] = {

    Transformation {

      (pt: Point[_3D]) =>
      {
        var meanPoint = EuclideanVector3D(0, 0, 0)
        var i = 0
        while (i < dataItems.size) {
          meanPoint += dataItems(i).transformation(pt).toVector
          i += 1
        }
        (meanPoint / dataItems.size).toPoint

      }
    }
  }
}




object DataCollectionOfMultiMesh {



  private def refdomain(refjoint: JointTriangleMesh[_3D]): UnstructuredPointsDomain[_3D] = {

    val refrot1 =refjoint.poseIds(0).map{pi=>
      val p=refjoint.objects(0).pointSet.point(pi)
      TranslationTransform(EuclideanVector3D(-0.01, 0.0, 0.0)).apply(p)}.toIndexedSeq


    val refrot2 =refjoint.poseIds(1).map{pi=>
      val p=refjoint.objects(1).pointSet.point(pi)
      TranslationTransform(EuclideanVector3D(-0.01, 0.0, 0.0)).apply(p)}.toIndexedSeq

    val refrot3 =refjoint.poseIds(2).map{pi=>
      val p=refjoint.objects(2).pointSet.point(pi)
      TranslationTransform(EuclideanVector3D(-0.01, 0.0, 0.0)).apply(p)}.toIndexedSeq

    val reftrans1 = IndexedSeq(TranslationTransform(EuclideanVector3D(0.0, 0.0, 0.01)).apply(refjoint.objects(0).pointSet.point(refjoint.poseIds(0)(0))))

    val reftrans2 = IndexedSeq(TranslationTransform(EuclideanVector3D(0.0, 0.0, 0.01)).apply(refjoint.objects(1).pointSet.point(refjoint.poseIds(1)(0))))

    val reftrans3 = IndexedSeq(TranslationTransform(EuclideanVector3D(0.0, 0.0, 0.01)).apply(refjoint.objects(2).pointSet.point(refjoint.poseIds(2)(0))))

    val d: IndexedSeq[Point[_3D]] = refrot1++ refrot2++ refrot3++ reftrans1++ reftrans2 ++ reftrans3++
      refjoint.objects(0).pointSet.points.toIndexedSeq ++
      refjoint.objects(1).pointSet.points.toIndexedSeq ++ refjoint.objects(2).pointSet.points.toIndexedSeq

    UnstructuredPointsDomain[_3D](d)

  }


  private def shapePosedmainFromJointTriangleMesh(refjoint: JointTriangleMesh[_3D], jointtrianglemesh: JointTriangleMesh[_3D]): UnstructuredPointsDomain[_3D] = {

    val refland1=refjoint.poseIds(0).map(pi=>Landmark[_3D](pi.id.toString,refjoint.objects(0).pointSet.point(pi)))
    val targland1=refjoint.poseIds(0).map(pi=>Landmark[_3D](pi.id.toString,jointtrianglemesh.objects(0).pointSet.point(pi)))
    val t1 = LandmarkRegistration.rigid3DLandmarkRegistration(targland1, refland1, refjoint.rotCenters(0))


    val rot1 =refjoint.poseIds(0).map{pi=>
      val p=refjoint.objects(0).pointSet.point(pi)
      t1.rotation.inverse.apply(p)}.toIndexedSeq

    val trans1 = IndexedSeq(Point3D(t1.translation.inverse.parameters(0),t1.translation.inverse.parameters(1),t1.translation.inverse.parameters(2)))

    val mesh1=jointtrianglemesh.objects(0).transform(t1)

    val refland2=refjoint.poseIds(1).map(pi=>Landmark[_3D](pi.id.toString,refjoint.objects(1).pointSet.point(pi)))
    val targland2=refjoint.poseIds(1).map(pi=>Landmark[_3D](pi.id.toString,jointtrianglemesh.objects(1).pointSet.point(pi)))
    val t2 = LandmarkRegistration.rigid3DLandmarkRegistration(targland2, refland2, refjoint.rotCenters(1))


    val rot2 =refjoint.poseIds(1).map{pi=>
      val p=refjoint.objects(1).pointSet.point(pi)
      t2.rotation.inverse.apply(p)}.toIndexedSeq

    val trans2 = IndexedSeq(Point3D(t2.translation.inverse.parameters(0),t1.translation.inverse.parameters(1),t1.translation.inverse.parameters(2)))

    val mesh2=jointtrianglemesh.objects(1).transform(t2)


    val refland3=refjoint.poseIds(2).map(pi=>Landmark[_3D](pi.id.toString,refjoint.objects(2).pointSet.point(pi)))
    val targland3=refjoint.poseIds(2).map(pi=>Landmark[_3D](pi.id.toString,jointtrianglemesh.objects(2).pointSet.point(pi)))
    val t3 = LandmarkRegistration.rigid3DLandmarkRegistration(targland3, refland3, refjoint.rotCenters(2))


    val rot3 =refjoint.poseIds(2).map{pi=>
      val p=refjoint.objects(2).pointSet.point(pi)
      t3.rotation.inverse.apply(p)}.toIndexedSeq

    val trans3 = IndexedSeq(Point3D(t3.translation.inverse.parameters(0),t1.translation.inverse.parameters(1),t1.translation.inverse.parameters(2)))

    val mesh3=jointtrianglemesh.objects(2).transform(t3)



    // MeshIO.writeMesh(refjoint.objects(0).transform(t1.rotation.inverse), new File("E:\\PhD folders\\Tetrahedral mesh\\test meshe sfor joinr models to deleted\\mesh1_" + mesh1.area + ".stl"))
    // MeshIO.writeMesh(refjoint.objects(1).transform(t2.rotation.inverse), new File("E:\\PhD folders\\Tetrahedral mesh\\test meshe sfor joinr models to deleted\\mesh2_" + mesh2.area + ".stl"))

    val d = rot1 ++ rot2++ rot3 ++ trans1++ trans2++ trans3++ mesh1.pointSet.points.toIndexedSeq++
      mesh2.pointSet.points.toIndexedSeq ++ mesh3.pointSet.points.toIndexedSeq

    UnstructuredPointsDomain(d)
  }


  def fromMeshSequence(referenceMesh: JointTriangleMesh[_3D], registeredMeshes: Seq[JointTriangleMesh[_3D]])(implicit rng: Random): (Option[DataCollectionOfMultiMesh], Seq[Throwable]) = {
    val ref = refdomain(referenceMesh)
    val reg = registeredMeshes.map { m => shapePosedmainFromJointTriangleMesh(referenceMesh, m) }

    val (transformations, errors) = DataUtils.partitionSuccAndFailedTries(reg.map(DataUtils.multimeshToTransformation(ref, _)))
    val dc = DataCollectionOfMultiMesh(referenceMesh, transformations.map(DataItem("from multi mesh", _)))
    if (dc.size > 0) (Some(dc), errors) else (None, errors)
  }

  /**
    * Builds a dc instance from a reference mesh volume and a directory containing meshe volumes in correspondence with the reference.
    * Only vtk and stl meshes are currently supported.
    *
    * @return a data collection containing the valid elements as well as the list of errors for invalid items.
    */
  def fromMeshDirectory(referenceMesh: JointTriangleMesh[_3D], meshDirectory1: File, meshDirectory2: File, meshDirectory3: File)(implicit rng: Random): (Option[DataCollectionOfMultiMesh], Seq[Throwable]) = {
    val meshFileNames1 = meshDirectory1.listFiles().toSeq.filter(fn => fn.getAbsolutePath.endsWith(".vtk") || fn.getAbsolutePath.endsWith(".stl"))
    val meshFileNames2 = meshDirectory2.listFiles().toSeq.filter(fn => fn.getAbsolutePath.endsWith(".vtk") || fn.getAbsolutePath.endsWith(".stl"))
    val meshFileNames3 = meshDirectory3.listFiles().toSeq.filter(fn => fn.getAbsolutePath.endsWith(".vtk") || fn.getAbsolutePath.endsWith(".stl"))

    val (meshes1, ioErrors1) = DataUtils.partitionSuccAndFailedTries(for (meshFn <- meshFileNames1) yield {
      MeshIO.readMesh(meshFn).map(m => TriangleMesh3D(m.pointSet, referenceMesh.objects.apply(0).triangulation))

    })

    val (meshes2, ioErrors2) = DataUtils.partitionSuccAndFailedTries(for (meshFn <- meshFileNames2) yield {
      MeshIO.readMesh(meshFn).map(m => TriangleMesh3D(m.pointSet, referenceMesh.objects.apply(1).triangulation))

    })

    val (meshes3, ioErrors3) = DataUtils.partitionSuccAndFailedTries(for (meshFn <- meshFileNames3) yield {
      MeshIO.readMesh(meshFn).map(m => TriangleMesh3D(m.pointSet, referenceMesh.objects.apply(2).triangulation))

    })

    val meshes = for (i <- 0 to meshes1.size - 1) yield {
      val d = UnstructuredPointsDomain(meshes1.apply(i).pointSet.points.toIndexedSeq ++ meshes2.apply(i).pointSet.points.toIndexedSeq ++ meshes3.apply(i).pointSet.points.toIndexedSeq)
      JointTriangleMesh3D(3, d, List(meshes1.apply(i), meshes1.apply(i), meshes1.apply(i)),referenceMesh.poseIds,referenceMesh.rotCenters)
    }

    val (dc, meshErrors) = fromMeshSequence(referenceMesh, meshes)
    (dc, ioErrors1 ++ ioErrors2 ++ ioErrors3 ++ meshErrors)
  }

}

