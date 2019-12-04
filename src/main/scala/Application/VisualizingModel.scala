package Application

import java.io.File

import DmoGpm.Model
import DmoGpm.io.StatisticalModelIO
import breeze.linalg.DenseVector
import scalismo.ui.api.ScalismoUI

object VisualizingModel {

  def main(args: Array[String]): Unit = {
    scalismo.initialize()
    val ui = ScalismoUI()



    val model = StatisticalModelIO.readStatisticalMultMeshModel(new File("data\\lollipop\\JointModel.h5")).get

    val l=List(-3.0,-1.0,1.0,3.0)

    sampling(0,model,l,ui)




  }

  def sampling(mode: Int, model: Model, coefs:List[Double], ui: ScalismoUI) {

    var v = DenseVector.zeros[Double](model.rank)

    for (c <- coefs) {
      v(mode) = c
      ui.show(model.instance(v).objects(0), +c + "std object 1")
      ui.show(model.instance(v).objects(1), +c + "std object 2")
      ui.show(model.instance(v).objects(2), +c + "std object 2")
    }
  }


}