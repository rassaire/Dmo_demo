package DmoGpm
import breeze.stats.distributions.Gaussian
import breeze.linalg.{DenseMatrix, DenseVector}
import scalismo.common._
import scalismo.geometry._
import scalismo.numerics.PivotedCholesky.StoppingCriterion
import DmoGpm.{JointTriangleMesh, JointTriangleMesh3D}
import scalismo.registration.{GaussianProcessTransformationSpace, LandmarkRegistration, RigidTransformation, TranslationTransform}
import scalismo.statisticalmodel.{DiscreteLowRankGaussianProcess, LowRankGaussianProcess}
import scalismo.statisticalmodel.DiscreteLowRankGaussianProcess.{Eigenpair => DiscreteEigenpair, _}
import scalismo.statisticalmodel.LowRankGaussianProcess.Eigenpair
import scalismo.numerics.PivotedCholesky.RelativeTolerance
import DmoGpm.dataset.DataCollectionOfMultiMesh
import scalismo.utils.Random

import scala.util.{Failure, Success, Try}

case class Model private  (referenceJointMesh: JointTriangleMesh[_3D], gp: DiscreteLowRankGaussianProcess[_3D, UnstructuredPointsDomain[_3D], EuclideanVector[_3D]]) {

  /** @see [[scalismo.statisticalmodel.DiscreteLowRankGaussianProcess.rank]] */

  val rank = gp.rank
  private val GPinterpolate = GaussianProcessTransformationSpace(gp.interpolate(NearestNeighborInterpolator()))
  /**
    * The mean shape
    * @see [[DiscreteLowRankGaussianProcess.mean]]
    */
  lazy val mean: JointTriangleMesh[_3D] = warpReferences(DenseVector.zeros[Double](rank))


  private val refrot1 =referenceJointMesh.poseIds(0).map{pi=>
    val p=referenceJointMesh.objects(0).pointSet.point(pi)
    TranslationTransform(EuclideanVector3D(-0.01, 0.0, 0.0)).apply(p)}.toIndexedSeq


  private  val refrot2 =referenceJointMesh.poseIds(1).map{pi=>
    val p=referenceJointMesh.objects(1).pointSet.point(pi)
    TranslationTransform(EuclideanVector3D(-0.01, 0.0, 0.0)).apply(p)}.toIndexedSeq

  private val refrot3 =referenceJointMesh.poseIds(2).map{pi=>
    val p=referenceJointMesh.objects(2).pointSet.point(pi)
    TranslationTransform(EuclideanVector3D(-0.01, 0.0, 0.0)).apply(p)}.toIndexedSeq

  private val reftrans1 = IndexedSeq(TranslationTransform(EuclideanVector3D(0.0, 0.0, 0.01)).apply(referenceJointMesh.objects(0).pointSet.point(referenceJointMesh.poseIds(0)(0))))

  private val reftrans2 = IndexedSeq(TranslationTransform(EuclideanVector3D(0.0, 0.0, 0.01)).apply(referenceJointMesh.objects(1).pointSet.point(referenceJointMesh.poseIds(1)(0))))

  private val reftrans3 = IndexedSeq(TranslationTransform(EuclideanVector3D(0.0, 0.0, 0.01)).apply(referenceJointMesh.objects(2).pointSet.point(referenceJointMesh.poseIds(2)(0))))


  val meanDomain: UnstructuredPointsDomain[_3D]=UnstructuredPointsDomain[_3D](refrot1++ refrot2++ refrot3++ reftrans1++ reftrans2 ++ reftrans3++
    referenceJointMesh.objects(0).pointSet.points.toIndexedSeq ++ referenceJointMesh.objects(1).pointSet.points.toIndexedSeq ++
    referenceJointMesh.objects(2).pointSet.points.toIndexedSeq).transform(GPinterpolate.transformForParameters(DenseVector.zeros[Double](rank)))
  /**
    * The covariance between two points of the  mesh with given point id.
    * @see [[DiscreteLowRankGaussianProcess.cov]]
    */
  def cov(ptId1: PointId, ptId2: PointId) = gp.cov(ptId1, ptId2)

  /**
    * draws a random shape.
    * @see [[DiscreteLowRankGaussianProcess.sample]]
    */
  def sample()(implicit rand: Random) = {
    val standardNormal = Gaussian(0, 1)(rand.breezeRandBasis)
    val coeffs = standardNormal.sample(rank)
    warpReferences(DenseVector(coeffs.toArray))
  }

  /**
    * returns the probability density for an instance of the model
    * @param instanceCoefficients coefficients of the instance in the model. For shapes in correspondence, these can be obtained using the coefficients method
    *
    */

  def pdf(instanceCoefficients: DenseVector[Double]): Double = {
    val disVecField = gp.instance(instanceCoefficients)
    gp.pdf(disVecField)
  }

  /**
    * returns a shape that corresponds to a linear combination of the basis functions with the given coefficients c.
    *  @see [[DiscreteLowRankGaussianProcess.instance]]
    */
  def instance(c: DenseVector[Double]): JointTriangleMesh[_3D] = warpReferences(c)

  /**
    *  Returns a marginal StatisticalMeshVolumeModel, modelling deformations only on the chosen points of the reference
    *
    *  This method proceeds by clipping the reference mesh to keep only the indicated point identifiers, and then marginalizing the
    *  GP over those points. Notice that when clipping, not all indicated point ids will be part of the clipped mesh volume, as some points may not belong
    *  to any cells anymore. Therefore 2 behaviours are supported by this method :
    *
    *  1- in case some of the indicated pointIds remain after clipping and do form a mesh, a marginal model is returned only for those points
    *  2- in case none of the indicated points remain (they are not meshed), a reference mesh voulme with all indicated point Ids and no cells is constructed and a marginal
    *  over this new reference is returned
    *
    * @see [[DiscreteLowRankGaussianProcess.marginal]]
    */

  /*def marginal(ptIds: IndexedSeq[PointId]) = {
    val clippedReference = referenceMeshVolume.operations.clip(p => { !ptIds.contains(referenceMesh.pointSet.findClosestPoint(p).id) })
    // not all of the ptIds remain in the reference after clipping, since their cells might disappear
    val remainingPtIds = clippedReference.pointSet.points.map(p => referenceMesh.pointSet.findClosestPoint(p).id).toIndexedSeq
    if (remainingPtIds.isEmpty) {
      val newRef = TriangleMesh3D(UnstructuredPointsDomain(ptIds.map(id => referenceMesh.pointSet.point(id)).toIndexedSeq), TriangleList(IndexedSeq[TriangleCell]()))
      val marginalGP = gp.marginal(ptIds.toIndexedSeq)
      StatisticalMeshModel(newRef, marginalGP)
    } else {
      val marginalGP = gp.marginal(remainingPtIds)
      StatisticalMeshModel(clippedReference, marginalGP)
    }
  }*/



  private def warpReferences(c: DenseVector[Double]): JointTriangleMesh[_3D] = {
    val registrationTransformation = GPinterpolate.transformForParameters(c)

    val currentmesh1 = referenceJointMesh.objects.apply(0).transform(registrationTransformation)
    val currentmesh2 = referenceJointMesh.objects.apply(1).transform(registrationTransformation)
    val currentmesh3 = referenceJointMesh.objects.apply(2).transform(registrationTransformation)

    val trans1p= reftrans1.map(m =>  registrationTransformation.apply(m))
    val trans1=TranslationTransform(trans1p(0).toVector)

    val rot1refland = refrot1.map(m => new Landmark[_3D](referenceJointMesh.objects.apply(0).pointSet.findClosestPoint(m).id.toString, registrationTransformation.apply(m)))
    val rot1tagland = refrot1.map(m => new Landmark[_3D](referenceJointMesh.objects.apply(0).pointSet.findClosestPoint(m).id.toString, m))
    val rot1 = LandmarkRegistration.rigid3DLandmarkRegistration(rot1tagland, rot1refland, referenceJointMesh.rotCenters(0))

    val mesh1 = currentmesh1.transform(RigidTransformation(rot1.rotation, trans1)) //(rot1.rotation).transform(trans1.translation)



    val trans2p= reftrans2.map(m =>  registrationTransformation.apply(m))
    val trans2=TranslationTransform(trans2p(0).toVector)

    val rot2refland = refrot2.map(m => new Landmark[_3D](referenceJointMesh.objects.apply(1).pointSet.findClosestPoint(m).id.toString, registrationTransformation.apply(m)))
    val rot2tagland = refrot2.map(m => new Landmark[_3D](referenceJointMesh.objects.apply(1).pointSet.findClosestPoint(m).id.toString, m))
    val rot2 = LandmarkRegistration.rigid3DLandmarkRegistration(rot2tagland, rot2refland, referenceJointMesh.rotCenters(1))

    val mesh2 = currentmesh2.transform(RigidTransformation(rot2.rotation, trans2)) //(rot1.rotation).transform(trans1.translation)






    val trans3p= reftrans3.map(m =>  registrationTransformation.apply(m))
    val trans3=TranslationTransform(trans3p(0).toVector)

    val rot3refland = refrot3.map(m => new Landmark[_3D](referenceJointMesh.objects.apply(2).pointSet.findClosestPoint(m).id.toString, registrationTransformation.apply(m)))
    val rot3tagland = refrot3.map(m => new Landmark[_3D](referenceJointMesh.objects.apply(2).pointSet.findClosestPoint(m).id.toString, m))
    val rot3 = LandmarkRegistration.rigid3DLandmarkRegistration(rot3tagland, rot3refland, referenceJointMesh.rotCenters(1))

    val mesh3 = currentmesh3.transform(RigidTransformation(rot3.rotation, trans3)) //(rot1.rotation).transform(trans1.translation)

    JointTriangleMesh3D(3, UnstructuredPointsDomain(mesh1.pointSet.points.toIndexedSeq ++ mesh2.pointSet.points.toIndexedSeq ++ mesh3.pointSet.points.toIndexedSeq),
      List(mesh1, mesh2, mesh3),referenceJointMesh.poseIds,referenceJointMesh.rotCenters)

  }

}

object Model {






  /**
    * creates a StatisticalMeshModel by discretizing the given Gaussian Process on the points of the reference mesh.
    */
  def apply(referenceMesh: JointTriangleMesh[_3D], gp: LowRankGaussianProcess[_3D, EuclideanVector[_3D]]): Model = {
    val discreteGp = DiscreteLowRankGaussianProcess(referenceMesh.pointSet, gp)
    new Model(referenceMesh, discreteGp)
  }

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

  /**
    * creates a StatisticalMeshModel from vector/matrix representation of the mean, variance and basis matrix.
    *
    * @see [[DiscreteLowRankGaussianProcess.apply(FiniteDiscreteDomain, DenseVector[Double], DenseVector[Double], DenseMatrix[Double]]
    */
  def apply(referenceMesh: JointTriangleMesh[_3D],
            meanVector: DenseVector[Double],
            variance: DenseVector[Double],
            basisMatrix: DenseMatrix[Double]) = {
    val gp =  DiscreteLowRankGaussianProcess[_3D, UnstructuredPointsDomain[_3D], EuclideanVector[_3D]](refdomain(referenceMesh), meanVector, variance, basisMatrix)
    new Model(referenceMesh, gp)
  }

  /**
    * Creates a new DiscreteLowRankGaussianProcess, where the mean and covariance matrix are estimated from the given transformations.
    *
    */
  def createUsingPCA(referenceMesh: JointTriangleMesh[_3D], fields: Seq[Field[_3D, EuclideanVector[_3D]]]): Model = {
    val c:StoppingCriterion = RelativeTolerance(0.0)
    val dgp:DiscreteLowRankGaussianProcess[_3D, UnstructuredPointsDomain[_3D], EuclideanVector[_3D]] = DiscreteLowRankGaussianProcess.createUsingPCA(refdomain(referenceMesh),fields,c)
    new Model(referenceMesh, dgp)
  }

  /**
    * Returns a PCA model with given reference mesh and a set of items in correspondence.
    * All points of the reference mesh are considered for computing the PCA
    */
  def createUsingPCA(dc: DataCollectionOfMultiMesh): Try[Model] = {
    if (dc.size < 3) return Failure(new Throwable(s"A data collection with at least 3 transformations is required to build a PCA Model (only ${dc.size} were provided)"))

    val fields = dc.dataItems.map { i =>
      Field[_3D, EuclideanVector[_3D]](i.transformation.domain, p => i.transformation(p) - p)
    }
    Success(createUsingPCA(dc.shapeReference, fields))
  }

}