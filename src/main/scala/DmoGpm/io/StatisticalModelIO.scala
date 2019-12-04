package DmoGpm.io

import java.io._
import java.util.Calendar

import breeze.linalg.{DenseMatrix, DenseVector}
import ncsa.hdf.`object`._
import scalismo.common.UnstructuredPointsDomain.Create.CreateUnstructuredPointsDomain3D
import scalismo.common.{PointId, UnstructuredPointsDomain}
import scalismo.geometry.{EuclideanVector3D, Landmark, Point, _3D}
import scalismo.io._
import scalismo.io.StatismoIO.StatismoModelType.StatismoModelType
import scalismo.mesh
import scalismo.mesh.TriangleMesh._
import scalismo.mesh.{TriangleCell, TriangleList, TriangleMesh, TriangleMesh3D}
import DmoGpm.{JointTriangleMesh, JointTriangleMesh3D}
import scalismo.registration.TranslationTransform
import scalismo.statisticalmodel.StatisticalMeshModel
import scalismo.mesh.{TetrahedralCell, TetrahedralList, TetrahedralMesh, TetrahedralMesh3D}
import DmoGpm.Model

import scala.util.{Failure, Success, Try}

object StatisticalModelIO {

  /**
    * Reads a statistical mesh model. The file type is determined
    * based on the extension. Currently on the Scalismo format (.h5)
    * is supported.
    *
    * @param file The statismo file
    * @return A StatisticalMeshModel or the Failure
    */


  /**
    * Reads a statistical mesh volume model. The file type is determined
    * based on the extension. Currently on the Scalismo format (.h5)
    * is supported.
    *
    * @param file The statismo file
    * @return A StatisticalMeshVolumeModel or the Failure
    */




  /**
    * Reads a statistical multi mesh model. The file type is determined
    * based on the extension. Currently on the Scalismo format (.h5)
    * is supported.
    *
    * @param file The statismo file
    * @return A Model or the Failure
    */
  def readStatisticalMultMeshModel(file: File): Try[Model] = {
    // currently, we support only the statismo format
    StatismoIO.readStatismoMultiMeshModel(file, "/")
  }

  /**
    * Writes a statistical mesh model. The file type is determined
    * based on the extension. Currently on the Scalismo format (.h5)
    * is supported.
    *
    * @param model The statistical model
    * @param file The file to which the model is written
    * @return In case of Failure, the Failure is returned.
    */

  /**
    * Writes a statistical mesh Volume model. The file type is determined
    * based on the extension. Currently on the Scalismo format (.h5)
    * is supported.
    *
    * @param model The statistical mesh volume model
    * @param file The file to which the model is written
    * @return In case of Failure, the Failure is returned.
    */


  /**
    * Writes a statistical multi meshmodel. The file type is determined
    * based on the extension. Currently on the Scalismo format (.h5)
    * is supported.
    *
    * @param model The statistical multi mesh model
    * @param file The file to which the model is written
    * @return In case of Failure, the Failure is returned.
    */
  def writeModel(model: Model, file: File): Try[Unit] = {
    // currently, we support only the statismo format
    StatismoIO.writeStatismoMultiMeshModel(model, file, "/")
  }
}

object StatismoIO {
  //
  //  object StatismoModelType extends Enumeration {
  //    type StatismoModelType = Value
  //    val Polygon_Mesh, Unknown = Value
  //
  //    def fromString(s: String): Value = {
  //      s match {
  //        case "POLYGON_MESH_MODEL" => Polygon_Mesh
  //        case _ => Unknown
  //      }
  //    }
  //  }
  //
  //  type ModelCatalog = Seq[CatalogEntry]
  //
  //  case class CatalogEntry(name: String, modelType: StatismoModelType, modelPath: String)
  //  object NoCatalogPresentException extends Exception
  //
  //  /**
  //    * List all models that are stored in the given hdf5 file.
  //    */
  //  def readModelCatalog(file: File): Try[ModelCatalog] = {
  //    import scala.collection.JavaConverters._
  //
  //    def flatten[A](xs: Seq[Try[A]]): Try[Seq[A]] = Try(xs.map(_.get))
  //
  //    for {
  //      h5file <- HDF5Utils.openFileForReading(file)
  //      catalogGroup <- if (h5file.exists("/catalog")) h5file.getGroup("/catalog") else Failure(NoCatalogPresentException)
  //      modelEntries = for (entryGroupObj <- catalogGroup.getMemberList.asScala.toSeq if entryGroupObj.isInstanceOf[Group]) yield {
  //        val entryGroup = entryGroupObj.asInstanceOf[Group]
  //        readCatalogEntry(h5file, entryGroup)
  //      }
  //      modelCatalog <- flatten(modelEntries)
  //    } yield {
  //      modelCatalog
  //    }
  //  }
  //
  //  private def readCatalogEntry(h5file: HDF5File, entryGroup: Group): Try[CatalogEntry] = {
  //    val name = entryGroup.getName
  //    for {
  //      location <- h5file.readString(entryGroup.getFullName + "/modelPath")
  //      modelType <- h5file.readString(entryGroup.getFullName + "/modelType")
  //    } yield {
  //      CatalogEntry(name, StatismoModelType.fromString(modelType), location)
  //    }
  //  }


  /**
    * Reads a statistical multi mesh model from a statismo file
    *
    * @param file      The statismo file
    * @param modelPath a path in the hdf5 file where the model is stored
    * @return
    */

  def readStatismoMultiMeshModel(file: File, modelPath: String = "/"): Try[Model] = {


    def extractOrthonormalPCABasisMatrix(pcaBasisMatrix: DenseMatrix[Double], pcaVarianceVector: DenseVector[Double]): DenseMatrix[Double] = {
      // this is an old statismo format, that has the pcaVariance directly stored in the PCA matrix,
      // i.e. pcaBasis = U * sqrt(lambda), where U is a matrix of eigenvectors and lambda the corresponding eigenvalues.
      // We recover U from it.

      val lambdaSqrt = pcaVarianceVector.map(l => math.sqrt(l))
      val lambdaSqrtInv = lambdaSqrt.map(l => if (l > 1e-8) 1.0f / l else 0f)

      // The following code is an efficient way to compute: pcaBasisMatrix * breeze.linalg.diag(lambdaSqrtInv)
      // (diag returns densematrix, so the direct computation would be very slow)
      val U = DenseMatrix.zeros[Double](pcaBasisMatrix.rows, pcaBasisMatrix.cols)
      for (i <- 0 until pcaBasisMatrix.cols) {
        U(::, i) := pcaBasisMatrix(::, i) * lambdaSqrtInv(i)
      }
      U
    }

    val modelOrFailure = for {
      h5file <- HDF5Utils.openFileForReading(file)

      representerName <- h5file.readStringAttribute(s"$modelPath/representer/", "name")
      // read mesh according to type given in representer
      mesh <- representerName match {
        case "multivtkPolyDataRepresenter" => readVTKMultiMeshFromRepresenterGroup(h5file, modelPath)
        case "multiitkMeshRepresenter" => readVTKMultiMeshFromRepresenterGroup(h5file, modelPath)
        case _ =>
          h5file.readStringAttribute(s"$modelPath/representer/", "datasetType") match {
            case Success("POLYGON_MESH") => readStandardMultiMeshFromRepresenterGroup(h5file, modelPath)
            case Success(datasetType) => Failure(new Exception(s"can only read model of datasetType POLYGON_MESH. Got $datasetType instead"))
            case Failure(t) => Failure(t)
          }
      }

      meanArray <- h5file.readNDArray[Float](s"$modelPath/model/mean")
      meanVector = DenseVector(meanArray.data.map(_.toDouble))
      pcaBasisArray <- h5file.readNDArray[Float](s"$modelPath/model/pcaBasis")
      majorVersion <- if (h5file.exists("/version/majorVersion")) h5file.readInt("/version/majorVersion")
      else {
        if (representerName == "multivtkPolyDataRepresenter" || representerName == "multiitkMeshRepresenter") Success(0)
        else Failure(new Throwable(s"no entry /version/majorVersion provided in statismo file."))
      }
      minorVersion <- if (h5file.exists("/version/minorVersion")) h5file.readInt("/version/minorVersion")
      else {
        if (representerName == "multivtkPolyDataRepresenter" || representerName == "multiitkMeshRepresenter") Success(8)
        else Failure(new Throwable(s"no entry /version/minorVersion provided in statismo file."))
      }
      pcaVarianceArray <- h5file.readNDArray[Float](s"$modelPath/model/pcaVariance")
      pcaVarianceVector = DenseVector(pcaVarianceArray.data.map(_.toDouble))
      pcaBasisMatrix = ndFloatArrayToDoubleMatrix(pcaBasisArray)
      pcaBasis <- (majorVersion, minorVersion) match {
        case (1, _) => Success(pcaBasisMatrix)
        case (0, 9) => Success(pcaBasisMatrix)
        case (0, 8) => Success(extractOrthonormalPCABasisMatrix(pcaBasisMatrix, pcaVarianceVector)) // an old statismo version
        case v => Failure(new Throwable(s"Unsupported version ${v._1}.${v._2}"))
      }

      _ <- Try {
        h5file.close()
      }
    } yield {
      // statismo stores the mean as the point position, not as a displacement on the reference.
      def flatten(v: IndexedSeq[Point[_3D]]) = DenseVector(v.flatten(pt => Array(pt(0), pt(1), pt(2))).toArray)

      val refpointsVec = flatten(mesh._1.points.toIndexedSeq)
      val meanDefVector = meanVector - refpointsVec

      Model(mesh._2, meanDefVector, pcaVarianceVector, pcaBasis)
    }

    modelOrFailure
  }


  object StatismoVersion extends Enumeration {
    type StatismoVersion = Value
    val v081, v090 = Value
  }

  import StatismoVersion._


  def writeStatismoMultiMeshModel(model: Model, file: File, modelPath: String = "/", statismoVersion: StatismoVersion = v090): Try[Unit] = {

    val discretizedMean = model.meanDomain.points.toIndexedSeq.flatten(_.toArray)
    val variance = model.gp.variance

    val pcaBasis = model.gp.basisMatrix.copy
    if (statismoVersion == v081) {
      // statismo 081 has the variance included in the pcaBasis
      for (i <- 0 until variance.length) {
        pcaBasis(::, i) *= math.sqrt(variance(i))
      }
    }
    val maybeError = for {
      h5file <- HDF5Utils.createFile(file)
      _ <- h5file.writeArray[Float](s"$modelPath/model/mean", discretizedMean.toArray.map(_.toFloat))
      _ <- h5file.writeArray[Float](s"$modelPath/model/noiseVariance", Array(0f))
      _ <- h5file.writeNDArray[Float](s"$modelPath/model/pcaBasis", NDArray(Array(pcaBasis.rows, pcaBasis.cols).map(_.toLong).toIndexedSeq, pcaBasis.t.flatten(false).toArray.map(_.toFloat)))
      _ <- h5file.writeArray[Float](s"$modelPath/model/pcaVariance", variance.toArray.map(_.toFloat))
      _ <- h5file.writeString(s"$modelPath/modelinfo/build-time", Calendar.getInstance.getTime.toString)
      group <- h5file.createGroup(s"$modelPath/representer")
      _ <- if (statismoVersion == v090) {
        for {
          _ <- writeRepresenterStatismov090_multiMesh(h5file, group, model, modelPath)
          _ <- h5file.writeInt("/version/majorVersion", 0)
          _ <- h5file.writeInt("/version/minorVersion", 9)
        } yield Success(())
      } else {
        for {
          _ <- writeRepresenterStatismov081_multiMesh(h5file, group, model, modelPath)
          _ <- h5file.writeInt("/version/majorVersion", 0)
          _ <- h5file.writeInt("/version/minorVersion", 8)
        } yield Success(())
      }
      _ <- h5file.writeString(s"$modelPath/modelinfo/modelBuilder-0/buildTime", Calendar.getInstance.getTime.toString)
      _ <- h5file.writeString(s"$modelPath/modelinfo/modelBuilder-0/builderName", "This is a useless info. The stkCore did not handle Model builder info at creation time.")
      _ <- h5file.createGroup(s"$modelPath/modelinfo/modelBuilder-0/parameters")
      _ <- h5file.createGroup(s"$modelPath/modelinfo/modelBuilder-0/dataInfo")
      _ <- Try {
        h5file.close()
      }
    } yield ()

    maybeError
  }


  private def writeRepresenterStatismov090_multiMesh(h5file: HDF5File, group: Group, model: Model, modelPath: String): Try[Unit] = {

    val cellArray1 = model.referenceJointMesh.objects(0).cells.map(_.ptId1.id) ++ model.referenceJointMesh.objects(0).cells.map(_.ptId2.id) ++ model.referenceJointMesh.objects(0).cells.map(_.ptId3.id)
    val pts1 = model.referenceJointMesh.objects(0).pointSet.points.toIndexedSeq.par.map(p => (p.toArray(0).toDouble, p.toArray(1).toDouble, p.toArray(2).toDouble))
    val pointArray1 = pts1.map(_._1) ++ pts1.map(_._2) ++ pts1.map(_._3)

    val cellArray2 = model.referenceJointMesh.objects(1).cells.map(_.ptId1.id) ++ model.referenceJointMesh.objects(1).cells.map(_.ptId2.id) ++ model.referenceJointMesh.objects(1).cells.map(_.ptId3.id)
    val pts2 = model.referenceJointMesh.objects(1).pointSet.points.toIndexedSeq.par.map(p => (p.toArray(0).toDouble, p.toArray(1).toDouble, p.toArray(2).toDouble))
    val pointArray2 = pts2.map(_._1) ++ pts2.map(_._2) ++ pts2.map(_._3)

    val cellArray3 = model.referenceJointMesh.objects(2).cells.map(_.ptId1.id) ++ model.referenceJointMesh.objects(2).cells.map(_.ptId2.id) ++ model.referenceJointMesh.objects(2).cells.map(_.ptId3.id)
    val pts3 = model.referenceJointMesh.objects(2).pointSet.points.toIndexedSeq.par.map(p => (p.toArray(0).toDouble, p.toArray(1).toDouble, p.toArray(2).toDouble))
    val pointArray3 = pts3.map(_._1) ++ pts3.map(_._2) ++ pts3.map(_._3)

    val cellArray1ids = model.referenceJointMesh.poseIds(0).map(pi => pi.id)
    val cellArray2ids = model.referenceJointMesh.poseIds(1).map(pi => pi.id)
    val cellArray3ids = model.referenceJointMesh.poseIds(2).map(pi => pi.id)
    val rotcenterspts = model.referenceJointMesh.rotCenters.toIndexedSeq.par.map(p => (p.toArray(0), p.toArray(1), p.toArray(2)))
    val rotcenters = rotcenterspts.map(_._1) ++ rotcenterspts.map(_._2) ++ rotcenterspts.map(_._3)


    for {
      _ <- h5file.writeStringAttribute(group.getFullName, "name", "multiitkStandardMeshRepresenter")
      _ <- h5file.writeStringAttribute(group.getFullName, "version/majorVersion", "0")
      _ <- h5file.writeStringAttribute(group.getFullName, "version/minorVersion", "9")
      _ <- h5file.writeStringAttribute(group.getFullName, "datasetType", "POLYGON_MESH")

      _ <- h5file.writeNDArray[Int](s"$modelPath/representer/cells1", NDArray(IndexedSeq(3, model.referenceJointMesh.objects(0).cells.size), cellArray1.toArray))
      _ <- h5file.writeNDArray[Float](s"$modelPath/representer/points1", NDArray(IndexedSeq(3, model.referenceJointMesh.objects(0).pointSet.points.size), pointArray1.toArray.map(_.toFloat)))

      _ <- h5file.writeNDArray[Int](s"$modelPath/representer/cells2", NDArray(IndexedSeq(3, model.referenceJointMesh.objects(1).cells.size), cellArray2.toArray))
      _ <- h5file.writeNDArray[Float](s"$modelPath/representer/points2", NDArray(IndexedSeq(3, model.referenceJointMesh.objects(1).pointSet.points.size), pointArray2.toArray.map(_.toFloat)))

      _ <- h5file.writeNDArray[Int](s"$modelPath/representer/cells3", NDArray(IndexedSeq(3, model.referenceJointMesh.objects(2).cells.size), cellArray3.toArray))
      _ <- h5file.writeNDArray[Float](s"$modelPath/representer/points3", NDArray(IndexedSeq(3, model.referenceJointMesh.objects(2).pointSet.points.size), pointArray3.toArray.map(_.toFloat)))

      _ <- h5file.writeNDArray[Int](s"$modelPath/representer/poseids1", NDArray(IndexedSeq(1, model.referenceJointMesh.poseIds(0).size), cellArray1ids.toArray))
      _ <- h5file.writeNDArray[Int](s"$modelPath/representer/poseids2", NDArray(IndexedSeq(1, model.referenceJointMesh.poseIds(1).size), cellArray2ids.toArray))
      _ <- h5file.writeNDArray[Int](s"$modelPath/representer/poseids3", NDArray(IndexedSeq(1, model.referenceJointMesh.poseIds(2).size), cellArray3ids.toArray))

      _ <- h5file.writeNDArray[Float](s"$modelPath/representer/rotcenters", NDArray(IndexedSeq(3, model.referenceJointMesh.rotCenters.size), rotcenters.toArray.map(_.toFloat)))
    } yield Success(())
  }


  private def writeRepresenterStatismov081_multiMesh(h5file: HDF5File, group: Group, model: Model, modelPath: String): Try[Unit] = {

    // we simply store the reference into a vtk file and store the file (the binary data) into the representer

    def refAsByteArray(ref: JointTriangleMesh[_3D]): Try[(Array[Byte], Array[Byte], Array[Byte], Array[Byte], Array[Byte], Array[Byte], Array[Byte])] = {
      val tmpfile1 = File.createTempFile("temp", ".vtk")
      tmpfile1.deleteOnExit()

      val tmpfile2 = File.createTempFile("temp", ".vtk")
      tmpfile2.deleteOnExit()

      val tmpfile3 = File.createTempFile("temp", ".vtk")
      tmpfile3.deleteOnExit()

      val tmpfile1id = File.createTempFile("temp", ".json")
      tmpfile1id.deleteOnExit()

      val tmpfile2id = File.createTempFile("temp", ".json")
      tmpfile2id.deleteOnExit()

      val tmpfile3id = File.createTempFile("temp", ".json")
      tmpfile3id.deleteOnExit()

      val tmpfilerotcenter = File.createTempFile("temp", ".json")
      tmpfilerotcenter.deleteOnExit()

      for {
        _ <- MeshIO.writeMesh(ref.objects(0), tmpfile1)
        _ <- MeshIO.writeMesh(ref.objects(1), tmpfile2)
        _ <- MeshIO.writeMesh(ref.objects(2), tmpfile3)
        _ <- LandmarkIO.writeLandmarksJson(ref.poseIds(0).map(pi => Landmark[_3D](pi.toString, ref.objects(0).pointSet.point(pi))), tmpfile1id)
        _ <- LandmarkIO.writeLandmarksJson(ref.poseIds(1).map(pi => Landmark[_3D](pi.toString, ref.objects(1).pointSet.point(pi))), tmpfile2id)
        _ <- LandmarkIO.writeLandmarksJson(ref.poseIds(2).map(pi => Landmark[_3D](pi.toString, ref.objects(2).pointSet.point(pi))), tmpfile3id)
        _ <- LandmarkIO.writeLandmarksJson(ref.rotCenters.map(p => Landmark[_3D](p.toVector.norm.toString, p)), tmpfile1id)
        ba <- readFileAsByteArray(tmpfile1, tmpfile2, tmpfile3, tmpfile1id, tmpfile2id, tmpfile3id, tmpfilerotcenter)
      } yield ba

    }

    def readFileAsByteArray(f1: File, f2: File, f3: File, f1id: File, f2id: File, f3id: File, frotcent: File): Try[(Array[Byte], Array[Byte], Array[Byte], Array[Byte], Array[Byte], Array[Byte], Array[Byte])] = {
      Try {
        val fileData1 = new Array[Byte](f1.length().toInt)
        val dis1 = new DataInputStream(new FileInputStream(f1))
        dis1.readFully(fileData1)
        dis1.close()

        val fileData2 = new Array[Byte](f2.length().toInt)
        val dis2 = new DataInputStream(new FileInputStream(f2))
        dis2.readFully(fileData2)
        dis2.close()

        val fileData3 = new Array[Byte](f3.length().toInt)
        val dis3 = new DataInputStream(new FileInputStream(f3))
        dis3.readFully(fileData3)
        dis3.close()

        val fileData1id = new Array[Byte](f1id.length().toInt)
        val dis1id = new DataInputStream(new FileInputStream(f1id))
        dis1id.readFully(fileData1id)
        dis1id.close()


        val fileData2id = new Array[Byte](f2id.length().toInt)
        val dis2id = new DataInputStream(new FileInputStream(f2id))
        dis2id.readFully(fileData2id)
        dis2id.close()

        val fileData3id = new Array[Byte](f3id.length().toInt)
        val dis3id = new DataInputStream(new FileInputStream(f3id))
        dis3id.readFully(fileData3id)
        dis3id.close()

        val fileDatarotcent = new Array[Byte](frotcent.length().toInt)
        val disrotcent = new DataInputStream(new FileInputStream(frotcent))
        disrotcent.readFully(fileDatarotcent)
        disrotcent.close()


        (fileData1, fileData2, fileData3, fileData1id, fileData2id, fileData3id, fileDatarotcent)
      }
    }

    for {
      _ <- h5file.writeStringAttribute(group.getFullName, "name", "multiitkMeshRepresenter")
      ba <- refAsByteArray(model.referenceJointMesh)
      _ <- h5file.writeNDArray[Byte](s"$modelPath/representer/reference1", NDArray(IndexedSeq(ba._1.length, 1), ba._1))
      _ <- h5file.writeNDArray[Byte](s"$modelPath/representer/reference2", NDArray(IndexedSeq(ba._2.length, 1), ba._2))
      _ <- h5file.writeNDArray[Byte](s"$modelPath/representer/reference3", NDArray(IndexedSeq(ba._3.length, 1), ba._3))
      _ <- h5file.writeNDArray[Byte](s"$modelPath/representer/poseids1", NDArray(IndexedSeq(ba._4.length, 1), ba._4))
      _ <- h5file.writeNDArray[Byte](s"$modelPath/representer/poseids2", NDArray(IndexedSeq(ba._5.length, 1), ba._5))
      _ <- h5file.writeNDArray[Byte](s"$modelPath/representer/poseids3", NDArray(IndexedSeq(ba._6.length, 1), ba._6))
      _ <- h5file.writeNDArray[Byte](s"$modelPath/representer/rotcenters", NDArray(IndexedSeq(ba._7.length, 1), ba._7))
    } yield Success(())
  }


  private def ndFloatArrayToDoubleMatrix(array: NDArray[Float])(implicit dummy: DummyImplicit, dummy2: DummyImplicit): DenseMatrix[Double] = {
    // the data in ndarray is stored row-major, but DenseMatrix stores it column major. We therefore
    // do switch dimensions and transpose
    DenseMatrix.create(array.dims(1).toInt, array.dims(0).toInt, array.data.map(_.toDouble)).t
  }

  private def ndIntArrayToIntMatrix(array: NDArray[Int]) = {
    // the data in ndarray is stored row-major, but DenseMatrix stores it column major. We therefore
    // do switch dimensions and transpose

    DenseMatrix.create(array.dims(1).toInt, array.dims(0).toInt, array.data).t
  }


  private def readStandardMultiMeshFromRepresenterGroup(h5file: HDF5File, modelPath: String):
  Try[(UnstructuredPointsDomain[_3D], JointTriangleMesh[_3D])] = {
    val mesh1 = for {
      vertArray1 <- h5file.readNDArray[Float](s"$modelPath/representer/points1").flatMap(vertArray =>
        if (vertArray.dims(0) != 3)
          Failure(new Exception("the representer points1 are not 3D points"))
        else
          Success(vertArray))
      vertMat = ndFloatArrayToDoubleMatrix(vertArray1)
      points = for (i <- 0 until vertMat.cols) yield Point(vertMat(0, i), vertMat(1, i), vertMat(2, i))
      cellArray <- h5file.readNDArray[Int](s"$modelPath/representer/cells1").flatMap(cellArray =>
        if (cellArray.dims(0) != 3)
          Failure(new Exception("the representer cells1 are not triangles"))
        else
          Success(cellArray))
      cellMat = ndIntArrayToIntMatrix(cellArray)
      cells = for (i <- 0 until cellMat.cols) yield TriangleCell(PointId(cellMat(0, i)), PointId(cellMat(1, i)), PointId(cellMat(2, i)))
      cellArray <- h5file.readNDArray[Int](s"$modelPath/representer/cells1")
    } yield TriangleMesh3D(UnstructuredPointsDomain(points), TriangleList(cells))

    val mesh2 = for {
      vertArray <- h5file.readNDArray[Float](s"$modelPath/representer/points2").flatMap(vertArray =>
        if (vertArray.dims(0) != 3)
          Failure(new Exception("the representer points2 are not 3D points"))
        else
          Success(vertArray))
      vertMat = ndFloatArrayToDoubleMatrix(vertArray)
      points = for (i <- 0 until vertMat.cols) yield Point(vertMat(0, i), vertMat(1, i), vertMat(2, i))
      cellArray <- h5file.readNDArray[Int](s"$modelPath/representer/cells2").flatMap(cellArray =>
        if (cellArray.dims(0) != 3)
          Failure(new Exception("the representer cells2 are not triangles"))
        else
          Success(cellArray))
      cellMat = ndIntArrayToIntMatrix(cellArray)
      cells = for (i <- 0 until cellMat.cols) yield TriangleCell(PointId(cellMat(0, i)), PointId(cellMat(1, i)), PointId(cellMat(2, i)))
      cellArray <- h5file.readNDArray[Int](s"$modelPath/representer/cells2")
    } yield TriangleMesh3D(UnstructuredPointsDomain(points), TriangleList(cells))

    val mesh3 = for {
      vertArray1 <- h5file.readNDArray[Float](s"$modelPath/representer/points3").flatMap(vertArray =>
        if (vertArray.dims(0) != 3)
          Failure(new Exception("the representer points3 are not 3D points"))
        else
          Success(vertArray))
      vertMat = ndFloatArrayToDoubleMatrix(vertArray1)
      points = for (i <- 0 until vertMat.cols) yield Point(vertMat(0, i), vertMat(1, i), vertMat(2, i))
      cellArray <- h5file.readNDArray[Int](s"$modelPath/representer/cells3").flatMap(cellArray =>
        if (cellArray.dims(0) != 3)
          Failure(new Exception("the representer cells3 are not triangles"))
        else
          Success(cellArray))
      cellMat = ndIntArrayToIntMatrix(cellArray)
      cells = for (i <- 0 until cellMat.cols) yield TriangleCell(PointId(cellMat(0, i)), PointId(cellMat(1, i)), PointId(cellMat(2, i)))
      cellArray <- h5file.readNDArray[Int](s"$modelPath/representer/cells3")
    } yield TriangleMesh3D(UnstructuredPointsDomain(points), TriangleList(cells))


    val poseids1 = for {
      cellArray <- h5file.readNDArray[Int](s"$modelPath/representer/poseids1").flatMap(cellArray =>
        if (cellArray.dims(0) != 1)
          Failure(new Exception("the representer ids are not pointid"))
        else
          Success(cellArray))
      poseids = for (i <- 0 until cellArray.data.size) yield PointId(cellArray.data(i))
    } yield poseids.toIndexedSeq


    val poseids2 = for {
      cellArray <- h5file.readNDArray[Int](s"$modelPath/representer/poseids2").flatMap(cellArray =>
        if (cellArray.dims(0) != 1)
          Failure(new Exception("the representer ids are not pointid"))
        else
          Success(cellArray))
      poseids = for (i <- 0 until cellArray.data.size) yield PointId(cellArray.data(i))
    } yield poseids.toIndexedSeq


    val poseids3 = for {
      cellArray <- h5file.readNDArray[Int](s"$modelPath/representer/poseids3").flatMap(cellArray =>
        if (cellArray.dims(0) != 1)
          Failure(new Exception("the representer ids are not pointid"))
        else
          Success(cellArray))
      poseids = for (i <- 0 until cellArray.data.size) yield PointId(cellArray.data(i))
    } yield poseids.toIndexedSeq


    val rotcenters = for {
      vertArray <- h5file.readNDArray[Float](s"$modelPath/representer/rotcenters").flatMap(vertArray =>
        if (vertArray.dims(0) != 3)
          Failure(new Exception("the representer points are not 3D points"))
        else
          Success(vertArray))
      vertMat = ndFloatArrayToDoubleMatrix(vertArray)
      points = for (i <- 0 until vertMat.cols) yield Point(vertMat(0, i), vertMat(1, i), vertMat(2, i))
    } yield points.toList


    val refrot1 = poseids1.get.map { pi =>
      val p = mesh1.get.pointSet.point(pi)
      TranslationTransform(EuclideanVector3D(-0.01, 0.0, 0.0)).apply(p)
    }.toIndexedSeq


    val refrot2 = poseids2.get.map { pi =>
      val p = mesh2.get.pointSet.point(pi)
      TranslationTransform(EuclideanVector3D(-0.01, 0.0, 0.0)).apply(p)
    }.toIndexedSeq

    val refrot3 = poseids3.get.map { pi =>
      val p = mesh3.get.pointSet.point(pi)
      TranslationTransform(EuclideanVector3D(-0.01, 0.0, 0.0)).apply(p)
    }.toIndexedSeq

    val reftrans1 = IndexedSeq(TranslationTransform(EuclideanVector3D(0.0, 0.0, 0.01)).apply(mesh1.get.pointSet.point(poseids1.get(0))))

    val reftrans2 = IndexedSeq(TranslationTransform(EuclideanVector3D(0.0, 0.0, 0.01)).apply(mesh2.get.pointSet.point(poseids2.get(0))))

    val reftrans3 = IndexedSeq(TranslationTransform(EuclideanVector3D(0.0, 0.0, 0.01)).apply(mesh3.get.pointSet.point(poseids3.get(0))))


    Try((UnstructuredPointsDomain[_3D](refrot1 ++ refrot2 ++ refrot3 ++ reftrans1 ++ reftrans2 ++ reftrans3 ++
      mesh1.get.pointSet.points.toIndexedSeq ++ mesh2.get.pointSet.points.toIndexedSeq ++
      mesh3.get.pointSet.points.toIndexedSeq), JointTriangleMesh3D(3, UnstructuredPointsDomain[_3D](mesh1.get.pointSet.points.toIndexedSeq ++ mesh2.get.pointSet.points.toIndexedSeq ++ mesh3.get.pointSet.points.toIndexedSeq),
      List(mesh1.get, mesh2.get, mesh3.get), List(poseids1.get, poseids2.get, poseids3.get), rotcenters.get)))
  }


  /*
   * reads the reference (a vtk file), which is stored as a byte array in the hdf5 file)
   */
  private def readVTKMeshFromRepresenterGroup(h5file: HDF5File, modelPath: String): Try[TriangleMesh[_3D]] = {
    for {
      rawdata <- h5file.readNDArray[Byte](s"$modelPath/representer/reference")
      vtkFile <- writeTmpFile(rawdata.data)
      triangleMesh <- MeshIO.readMesh(vtkFile)
    } yield triangleMesh
  }

  /*
 * reads the reference (a vtk file), which is stored as a byte array in the hdf5 file)
 */
  private def readVTKMeshVolumeFromRepresenterGroup(h5file: HDF5File, modelPath: String): Try[TetrahedralMesh[_3D]] = {
    for {
      rawdata <- h5file.readNDArray[Byte](s"$modelPath/representer/reference")
      vtkFile <- writeTmpFile(rawdata.data)
      tetrahedralMesh <- MeshIO.readTetrahedralMesh(vtkFile)
    } yield tetrahedralMesh
  }


  /*
* reads the joint reference (a vtk file), which is stored as a byte array in the hdf5 file)
*/
  private def readVTKMultiMeshFromRepresenterGroup(h5file: HDF5File, modelPath: String):
  Try[(UnstructuredPointsDomain[_3D], JointTriangleMesh[_3D])] = {
    val mesh1 = for {
      rawdata <- h5file.readNDArray[Byte](s"$modelPath/representer/reference1")
      vtkFile <- writeTmpFile(rawdata.data)
      triangleMesh <- MeshIO.readMesh(vtkFile)
    } yield triangleMesh

    val mesh2 = for {
      rawdata <- h5file.readNDArray[Byte](s"$modelPath/representer/reference2")
      vtkFile <- writeTmpFile(rawdata.data)
      triangleMesh <- MeshIO.readMesh(vtkFile)
    } yield triangleMesh

    val mesh3 = for {
      rawdata <- h5file.readNDArray[Byte](s"$modelPath/representer/reference3")
      vtkFile <- writeTmpFile(rawdata.data)
      triangleMesh <- MeshIO.readMesh(vtkFile)
    } yield triangleMesh


    val poseids1 = for {
      rawdata <- h5file.readNDArray[Byte](s"$modelPath/representer/poseids1")
      vtkFile <- writejsonTmpFile(rawdata.data)
      landmark <- LandmarkIO.readLandmarksJson[_3D](vtkFile)
    } yield landmark.map(l => PointId(l.id.toInt)).toIndexedSeq

    val poseids2 = for {
      rawdata <- h5file.readNDArray[Byte](s"$modelPath/representer/poseids2")
      vtkFile <- writejsonTmpFile(rawdata.data)
      landmark <- LandmarkIO.readLandmarksJson[_3D](vtkFile)
    } yield landmark.map(l => PointId(l.id.toInt)).toIndexedSeq

    val poseids3 = for {
      rawdata <- h5file.readNDArray[Byte](s"$modelPath/representer/poseids3")
      vtkFile <- writejsonTmpFile(rawdata.data)
      landmark <- LandmarkIO.readLandmarksJson[_3D](vtkFile)
    } yield landmark.map(l => PointId(l.id.toInt)).toIndexedSeq


    val rotcenters = for {
      rawdata <- h5file.readNDArray[Byte](s"$modelPath/representer/rotcenters")
      vtkFile <- writejsonTmpFile(rawdata.data)
      landmark <- LandmarkIO.readLandmarksJson[_3D](vtkFile)
    } yield landmark.map(l => l.point).toList


    val refrot1 = poseids1.get.map { pi =>
      val p = mesh1.get.pointSet.point(pi)
      TranslationTransform(EuclideanVector3D(-0.01, 0.0, 0.0)).apply(p)
    }.toIndexedSeq


    val refrot2 = poseids2.get.map { pi =>
      val p = mesh2.get.pointSet.point(pi)
      TranslationTransform(EuclideanVector3D(-0.01, 0.0, 0.0)).apply(p)
    }.toIndexedSeq

    val refrot3 = poseids3.get.map { pi =>
      val p = mesh3.get.pointSet.point(pi)
      TranslationTransform(EuclideanVector3D(-0.01, 0.0, 0.0)).apply(p)
    }.toIndexedSeq

    val reftrans1 = IndexedSeq(TranslationTransform(EuclideanVector3D(0.0, 0.0, 0.01)).apply(mesh1.get.pointSet.point(poseids1.get(0))))

    val reftrans2 = IndexedSeq(TranslationTransform(EuclideanVector3D(0.0, 0.0, 0.01)).apply(mesh2.get.pointSet.point(poseids2.get(0))))

    val reftrans3 = IndexedSeq(TranslationTransform(EuclideanVector3D(0.0, 0.0, 0.01)).apply(mesh3.get.pointSet.point(poseids3.get(0))))


    Try((UnstructuredPointsDomain[_3D](refrot1 ++ refrot2 ++ refrot3 ++ reftrans1 ++ reftrans2 ++ reftrans3 ++
      mesh1.get.pointSet.points.toIndexedSeq ++ mesh2.get.pointSet.points.toIndexedSeq ++
      mesh3.get.pointSet.points.toIndexedSeq), JointTriangleMesh3D(3, UnstructuredPointsDomain[_3D](mesh1.get.pointSet.points.toIndexedSeq ++ mesh2.get.pointSet.points.toIndexedSeq ++ mesh3.get.pointSet.points.toIndexedSeq),
      List(mesh1.get, mesh2.get, mesh3.get), List(poseids1.get, poseids2.get, poseids3.get), rotcenters.get)))
  }


  private def writeTmpFile(data: Array[Byte]): Try[File] = {
    val tmpfile = File.createTempFile("temp", ".vtk")
    tmpfile.deleteOnExit()

    Try {
      val stream = new DataOutputStream(new FileOutputStream(tmpfile))
      stream.write(data)
      stream.close()
    } map (_ => tmpfile)
  }


  private def writejsonTmpFile(data: Array[Byte]): Try[File] = {
    val tmpfile = File.createTempFile("temp", ".json")
    tmpfile.deleteOnExit()
    Try {
      val stream = new DataOutputStream(new FileOutputStream(tmpfile))
      stream.write(data)
      stream.close()
    } map (_ => tmpfile)
  }

  private def refdomain(refjoint: JointTriangleMesh[_3D]): UnstructuredPointsDomain[_3D] = {

    val refrot1 = refjoint.poseIds(0).map { pi =>
      val p = refjoint.objects(0).pointSet.point(pi)
      TranslationTransform(EuclideanVector3D(-0.01, 0.0, 0.0)).apply(p)
    }.toIndexedSeq


    val refrot2 = refjoint.poseIds(1).map { pi =>
      val p = refjoint.objects(1).pointSet.point(pi)
      TranslationTransform(EuclideanVector3D(-0.01, 0.0, 0.0)).apply(p)
    }.toIndexedSeq

    val refrot3 = refjoint.poseIds(2).map { pi =>
      val p = refjoint.objects(2).pointSet.point(pi)
      TranslationTransform(EuclideanVector3D(-0.01, 0.0, 0.0)).apply(p)
    }.toIndexedSeq

    val reftrans1 = IndexedSeq(TranslationTransform(EuclideanVector3D(0.0, 0.0, 0.01)).apply(refjoint.objects(0).pointSet.point(refjoint.poseIds(0)(0))))

    val reftrans2 = IndexedSeq(TranslationTransform(EuclideanVector3D(0.0, 0.0, 0.01)).apply(refjoint.objects(1).pointSet.point(refjoint.poseIds(1)(0))))

    val reftrans3 = IndexedSeq(TranslationTransform(EuclideanVector3D(0.0, 0.0, 0.01)).apply(refjoint.objects(2).pointSet.point(refjoint.poseIds(2)(0))))

    val d: IndexedSeq[Point[_3D]] = refrot1 ++ refrot2 ++ refrot3 ++ reftrans1 ++ reftrans2 ++ reftrans3 ++
      refjoint.objects(0).pointSet.points.toIndexedSeq ++
      refjoint.objects(1).pointSet.points.toIndexedSeq ++ refjoint.objects(2).pointSet.points.toIndexedSeq

    UnstructuredPointsDomain[_3D](d)

  }
}
