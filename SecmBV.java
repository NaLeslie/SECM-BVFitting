import com.comsol.model.*;
import com.comsol.model.util.*;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Paths;
import java.util.Calendar;
import java.util.LinkedList;
import java.util.Scanner;

/**
 *
 * @author Nathaniel Leslie
 */
public class SecmBV {
    static final double C_OX_BULK = 0.507;//mM
	static final double C_RED_BULK = 0.482;//mM
	static final double D_OX = 6.91E-10;//m^2/s
	static final double D_RED = 7.76E-10;//m^2/s
	static final double RG = 3.260;
	static final double UME_DIAM = 10;//um
	/**
     * nF/RT = 38.92 V<sup>-1</sup>
     */
    static final double F = 1.0 * 96485.3 / 8.314 / 294.0;
	/**
     * Transfer coefficient of the electrochemical system
     */
    static final double ALPHA = 0.58;//!!! This value has changed
	
	    /**
     * The values for the base-ten logarithm of the k-parameter that are to be tested to find the most suitable initial guess.
     */
    static final double[] TEST_LOG_K = new double[]{-7, -6.5, -6, -5.5, -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0};
    
    /**
     * The maximum number of iterations to be taken by the fitting method.
     */
    static final int MAX_ITERATIONS = 15;
    
    /**
     * The maximum value of lambda that can be obtained before giving up.
     */
    static final double MAX_LAMBDA = 10000.0;
    
    /**
     * The number of decimal points to which the L-parameter is to be considered.
     */
    static final int L_DECIMALS = 2;//in a
    
    /**
     * The perturbation that is made to the L-parameter
     */
    static final double L_PERTURB = 0.05;//in a
    
    /**
     * The number of decimal points to which the L-parameter is to be considered.
     */
    static final int ADJ_DECIMALS = 2;//in a
    
    /**
     * The perturbation that is made to the xa-parameter
     */
    static final double ADJ_X_PERTURB = 0.5;//in a
	
	/**
     * The perturbation that is made to the ya-parameter
     */
    static final double ADJ_Y_PERTURB = 0.5;//in a
    
    /**
     * The number of decimal points to which the overpotential is to be considered.
     */
    static final int E_DECIMALS = 3;//[E] = [V]
    
    /**
     * The perturbation that is to be introduced to the overpotential.
     */
    static final double E_PERTURB = 0.005;//[E] = [V]
    
	/**
     * The number of decimal points to which the base-ten logarithm of k is to be considered.
     */
    static final int LOGK_DECIMALS = 3;//[k] = [m/s]
    
    /**
     * The perturbation that is to be introduced to the base-ten logarithm of k.
     */
    static final double LOGK_PERTURB = 0.005;//[k] = [m/s]
	
    /**
     * The number of decimal points to which the R-parameter is to be considered.
     */
    static final int R_DECIMALS = 1;//in grid pixels
    
    /**
     * The perturbation that is made to the R-parameter
     */
    static final double R_PERTURB = 1.5;//in grid pixels
	
    /**
    * Stand-in for auto-generated method
    * @param reactivity_data_filename
    * @param l The L parameter
     * @param xadj The x-position adjustment
     * @param yadj The y-parameter adjustment
    * @param xcoords The x coordinates of each sampled point
    * @param ycoords The y coordinates of each sampled point
    * @return 
    */
	public static Model run(String reactivity_data_filename, double l, double xadj, double yadj, double[] xcoords, double[] ycoords, double[] redoxdirs) {
        Model model = ModelUtil.create("Model");
	

    //<MOD>
	String xpositions = toString(xcoords);
	String ypositions = toString(ycoords);
	String redoxdirections = toString(redoxdirs);
	//System.out.println(xpositions);
	//System.out.println(ypositions);
	//System.out.println(redoxdirections);

	String cwd = getCWD();

    model.modelPath(cwd);
	//</MOD>

    model.label("SECM_BV_img.mph");


    model.param().set("a", "UME_Diam / 2", "UME RADIUS");
    model.param().set("Rg", RG + "*a", "Normalized Glass RADIUS");
    model.param().set("UME_Diam", UME_DIAM + "[um]", "UME DIAMETER");
    model.param().set("X_UME", "Max_XY/2", "X-position of UME");
    model.param().set("Y_UME", "Max_XY/2", "Y-position of UME");
    //<MOD>
    model.param().set("Z_UME", l + "*a", "Z-position of UME");
    model.param().set("X_ADJ", xadj + "*a", "X-adjustment from the fit");
    model.param().set("Y_ADJ", yadj + "*a", "Y-adjustment from the fit");//</MOD>
    model.param().set("Max_XY", "a*350", "Maximum x,y bound");
    model.param().set("Max_Z", "Max_XY *1.0", "Maximum z bound");
    model.param().set("L_UME", "Max_Z - Z_UME", "Length of UME");
    model.param().set("Cox_bulk", C_OX_BULK + "[mM]", "Bulk concentration of oxidized species");
    model.param().set("Cred_bulk", C_RED_BULK + "[mM]", "Bulk concentration of reduced species");
    model.param().set("Dox", D_OX + "[m^2/s]", "Diffusion coefficient - oxidised");
	model.param().set("Dred", D_RED + "[m^2/s]", "Diffusion coefficient - reduced");
    model.param().set("k_UME", "1 [m/s]", "Rate constant at UME");
    model.param().set("k_PS", "10^logk [m/s]", "Rate constant at point source");
    model.param()
         .set("MaxMeshRadial_frac", "0.04", "Fraction of the radial features that will be the maximum mesh size");
    model.param().set("sl", "201", "grid side length");
    model.param().set("gridres", "0.1*a", "Grid resolution");
    model.param().set("logk", "-4.005", "base-10 logarithm of k/1m/s");
    model.param().set("xoffs", "1.351250E-05[m]");
    model.param().set("yoffs", "4.359250E-05[m]");
    model.param().set("factor", "10^(logk)/(1E-4)");
    model.param().set("o0r1", "0", "parameter if zero ME oxidizes, if 1 ME reduces");

    model.component().create("comp1", true);

    model.component("comp1").geom().create("geom1", 3);

    model.component("comp1").curvedInterior(false);

    model.result().table().create("tbl1", "Table");

    model.component("comp1").func().create("int1", "Interpolation");
    model.component("comp1").func().create("an1", "Analytic");
    model.component("comp1").func().create("an2", "Analytic");
    model.component("comp1").func("int1").set("source", "file");
    model.component("comp1").func("int1").label("rate");
    model.component("comp1").func("int1").set("funcs", new String[][]{{"red_rate_map", "1"}, {"ox_rate_map", "2"}});
    //<MOD>
    model.component("comp1").func("int1").set("filename", cwd + "/" + reactivity_data_filename);//</MOD>
    model.component("comp1").func("int1").set("nargs", 2);
    model.component("comp1").func("int1").set("interp", "neighbor");
    model.component("comp1").func("int1").set("extrap", "value");
    model.component("comp1").func("int1").set("fununit", new String[]{"m/s", "m/s"});
    model.component("comp1").func("int1").set("argunit", new String[]{"m", "m"});
    model.component("comp1").func("an1").label("ShiftedKred");
    model.component("comp1").func("an1").set("funcname", "kred_func");
    model.component("comp1").func("an1")
         .set("expr", "red_rate_map(x+xoffs-0.5*Max_XY+X_ADJ,y+yoffs-0.5*Max_XY+Y_ADJ)");
    model.component("comp1").func("an1").set("args", new String[]{"x", "y"});
    model.component("comp1").func("an1").set("fununit", "m/s");
    model.component("comp1").func("an1").set("argunit", new String[]{"m", "m"});
    model.component("comp1").func("an1")
         .set("plotargs", new String[][]{{"x", "Max_XY*0.4", "Max_XY*0.6"}, {"y", "Max_XY*0.4", "Max_XY*0.6"}});
    model.component("comp1").func("an2").label("ShiftedKox");
    model.component("comp1").func("an2").set("funcname", "kox_func");
    model.component("comp1").func("an2")
         .set("expr", "ox_rate_map(x+xoffs-0.5*Max_XY+X_ADJ,y+yoffs-0.5*Max_XY+Y_ADJ)");
    model.component("comp1").func("an2").set("args", new String[]{"x", "y"});
    model.component("comp1").func("an2").set("fununit", "m/s");
    model.component("comp1").func("an2").set("argunit", new String[]{"m", "m"});
    model.component("comp1").func("an2")
         .set("plotargs", new String[][]{{"x", "Max_XY*0.4", "Max_XY*0.6"}, {"y", "Max_XY*0.4", "Max_XY*0.6"}});

    model.component("comp1").mesh().create("mesh1");

    model.component("comp1").geom("geom1").label("Simulation cell");
    model.component("comp1").geom("geom1").create("blk1", "Block");
    model.component("comp1").geom("geom1").feature("blk1").label("Bath");
    model.component("comp1").geom("geom1").feature("blk1").set("size", new String[]{"Max_XY", "Max_XY", "Max_Z"});
    model.component("comp1").geom("geom1").create("cyl1", "Cylinder");
    model.component("comp1").geom("geom1").feature("cyl1").label("UME");
    model.component("comp1").geom("geom1").feature("cyl1").set("pos", new String[]{"X_UME", "Y_UME", "Z_UME"});
    model.component("comp1").geom("geom1").feature("cyl1").set("r", "a");
    model.component("comp1").geom("geom1").feature("cyl1").set("h", "L_UME");
    model.component("comp1").geom("geom1").create("cyl3", "Cylinder");
    model.component("comp1").geom("geom1").feature("cyl3").label("UME_Ins");
    model.component("comp1").geom("geom1").feature("cyl3").set("pos", new String[]{"X_UME", "Y_UME", "Z_UME"});
    model.component("comp1").geom("geom1").feature("cyl3").set("r", "Rg");
    model.component("comp1").geom("geom1").feature("cyl3").set("h", "L_UME");
    model.component("comp1").geom("geom1").create("blk2", "Block");
    model.component("comp1").geom("geom1").feature("blk2").active(false);
    model.component("comp1").geom("geom1").feature("blk2").label("box00");
    model.component("comp1").geom("geom1").feature("blk2").set("pos", new String[]{"42.25*a", "43.75*a", "-1E-5"});
    model.component("comp1").geom("geom1").feature("blk2").set("size", new String[]{"0.5*a", "0.5*a", "1E-5"});
    model.component("comp1").geom("geom1").run();
    model.component("comp1").geom("geom1").run("fin");

    model.view().create("view2", 2);
    model.view().create("view3", 2);

    model.component("comp1").physics().create("tds", "DilutedSpecies", "geom1");
    model.component("comp1").physics("tds").field("concentration").field("Cox");
    model.component("comp1").physics("tds").field("concentration").component(new String[]{"Cox", "Cred"});
    model.component("comp1").physics("tds").selection().set(1);
    model.component("comp1").physics("tds").create("fl1", "FluxBoundary", 2);
    model.component("comp1").physics("tds").feature("fl1").selection().set(12);
    model.component("comp1").physics("tds").create("fl2", "FluxBoundary", 2);
    model.component("comp1").physics("tds").feature("fl2").selection().set(3);
    model.component("comp1").physics("tds").create("conc1", "Concentration", 2);
    model.component("comp1").physics("tds").feature("conc1").selection().set(1, 2, 4, 5, 18);

    model.component("comp1").mesh("mesh1").create("ftet1", "FreeTet");
    model.component("comp1").mesh("mesh1").feature("ftet1").selection().geom("geom1", 3);
    model.component("comp1").mesh("mesh1").feature("ftet1").selection().set(1);
    model.component("comp1").mesh("mesh1").feature("ftet1").create("size1", "Size");
    model.component("comp1").mesh("mesh1").feature("ftet1").create("size2", "Size");
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size1").selection().geom("geom1", 2);
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size1").selection().set(12);
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size2").selection().geom("geom1", 2);
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size2").selection().set(3);

    model.component("comp1").probe().create("bnd1", "Boundary");
    model.component("comp1").probe("bnd1").selection().set(12);

    model.result().table("tbl1").label("Probe Table 1");

    model.component("comp1").view("view1").set("transparency", true);
    model.view("view2").axis().set("xmin", 397.5);
    model.view("view2").axis().set("xmax", 865);
    model.view("view2").axis().set("ymin", 494.375);
    model.view("view2").axis().set("ymax", 755.625);
    model.view("view2").axis().set("viewscaletype", "automatic");
    model.view("view3").axis().set("xmin", 41.69999694824219);
    model.view("view3").axis().set("xmax", 59.30000305175781);
    model.view("view3").axis().set("ymin", 42.8861083984375);
    model.view("view3").axis().set("ymax", 57.1138916015625);
    model.view("view3").axis().set("viewscaletype", "automatic");

    model.component("comp1").physics("tds").prop("EquationForm").set("form", "Stationary");
    model.component("comp1").physics("tds").prop("TransportMechanism").set("Convection", false);
    model.component("comp1").physics("tds").feature("cdm1")
         .set("D_Cox", new String[][]{{"Dox"}, {"0"}, {"0"}, {"0"}, {"Dox"}, {"0"}, {"0"}, {"0"}, {"Dox"}});
    model.component("comp1").physics("tds").feature("cdm1")
         .set("D_Cred", new String[][]{{"Dred"}, {"0"}, {"0"}, {"0"}, {"Dred"}, {"0"}, {"0"}, {"0"}, {"Dred"}});
    model.component("comp1").physics("tds").feature("init1")
         .set("initc", new String[][]{{"Cox_bulk"}, {"Cred_bulk"}});
    model.component("comp1").physics("tds").feature("fl1").set("species", new int[][]{{1}, {1}});
    model.component("comp1").physics("tds").feature("fl1")
         .set("J0", new String[][]{{"-k_UME*Cred*(o0r1-1) - k_UME*Cox*(o0r1)"}, {"k_UME*Cox*(o0r1) + k_UME*Cred*(o0r1-1)"}});
    model.component("comp1").physics("tds").feature("fl1").label("Flux_UME");
    model.component("comp1").physics("tds").feature("fl2").set("species", new int[][]{{1}, {1}});
    model.component("comp1").physics("tds").feature("fl2")
         .set("J0", new String[][]{{"(kox_func(x,y) * Cred - kred_func(x,y) * Cox)"}, {"(kred_func(x,y) * Cox - kox_func(x,y) * Cred)"}});
    model.component("comp1").physics("tds").feature("fl2").label("Flux_PS");
    model.component("comp1").physics("tds").feature("conc1").set("species", new int[][]{{1}, {1}});
    model.component("comp1").physics("tds").feature("conc1").set("c0", new String[][]{{"Cox_bulk"}, {"Cred_bulk"}});
    model.component("comp1").physics("tds").feature("conc1").label("Bulk");

    model.component("comp1").mesh("mesh1").feature("size").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("size").set("hmax", "80[um]");
    model.component("comp1").mesh("mesh1").feature("size").set("hmin", "10[um]");
    model.component("comp1").mesh("mesh1").feature("size").set("hgrad", 1.15);
    model.component("comp1").mesh("mesh1").feature("ftet1").set("optcurved", false);
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size1").label("UME");
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size1").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size1")
         .set("hmax", "a*MaxMeshRadial_frac*0.18");
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size1").set("hmaxactive", true);
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size1").set("hmin", 1.8E-5);
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size1").set("hminactive", false);
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size2").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size2").set("hgrad", "1.10");
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size2").set("hgradactive", true);
    model.component("comp1").mesh("mesh1").run();

    model.component("comp1").probe("bnd1").label("UME_FLUX");
    model.component("comp1").probe("bnd1").set("type", "integral");
    model.component("comp1").probe("bnd1").set("probename", "I_UME");
    model.component("comp1").probe("bnd1").set("expr", "tds.ntflux_Cred*F_const");
    model.component("comp1").probe("bnd1").set("unit", "A");
    model.component("comp1").probe("bnd1").set("descractive", true);
    model.component("comp1").probe("bnd1").set("descr", "Current at UME");
    model.component("comp1").probe("bnd1").set("table", "tbl1");
    model.component("comp1").probe("bnd1").set("window", "window1");

    model.study().create("std1");
    model.study("std1").create("param", "Parametric");
    model.study("std1").create("stat", "Stationary");

    model.sol().create("sol1");
    model.sol("sol1").study("std1");
    model.sol("sol1").attach("std1");
    model.sol("sol1").create("st1", "StudyStep");
    model.sol("sol1").create("v1", "Variables");
    model.sol("sol1").create("s1", "Stationary");
    model.sol("sol1").feature("s1").create("p1", "Parametric");
    model.sol("sol1").feature("s1").create("fc1", "FullyCoupled");
    model.sol("sol1").feature("s1").create("i1", "Iterative");
    model.sol("sol1").feature("s1").create("d1", "Direct");
    model.sol("sol1").feature("s1").feature("i1").create("mg1", "Multigrid");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").create("sl1", "SORLine");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").create("sl1", "SORLine");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").create("d1", "Direct");
    model.sol("sol1").feature("s1").feature().remove("fcDef");

    model.result().dataset().create("dset2", "Solution");
    model.result().dataset().create("int1", "Integral");
    model.result().dataset("dset2").set("probetag", "bnd1");
    model.result().dataset("int1").set("probetag", "bnd1");
    model.result().dataset("int1").set("data", "dset2");
    model.result().dataset("int1").selection().geom("geom1", 2);
    model.result().dataset("int1").selection().set(12);
    model.result().numerical().create("pev1", "EvalPoint");
    model.result().numerical("pev1").set("probetag", "bnd1");
    model.result().create("pg6", "PlotGroup2D");
    model.result().create("pg7", "PlotGroup2D");
    model.result().create("pg8", "PlotGroup3D");
    model.result().create("pg9", "PlotGroup3D");
    model.result().create("pg10", "PlotGroup3D");
    model.result().create("pg11", "PlotGroup3D");
    model.result().create("pg3", "PlotGroup1D");
    model.result("pg6").create("tbls1", "TableSurface");
    model.result("pg6").create("surf1", "Surface");
    model.result("pg6").feature("surf1").set("expr", "comp1.Cox");
    model.result("pg7").create("tbls1", "TableSurface");
    model.result("pg8").create("str1", "Streamline");
    model.result("pg8").feature("str1").create("col", "Color");
    model.result("pg9").create("surf1", "Surface");
    model.result("pg10").create("str1", "Streamline");
    model.result("pg10").feature("str1").create("col", "Color");
    model.result("pg10").feature("str1").feature("col").set("expr", "Cred");
    model.result("pg11").create("surf1", "Surface");
    model.result("pg11").feature("surf1").set("expr", "Cred");
    model.result("pg3").set("probetag", "window1");
    model.result("pg3").create("tblp1", "Table");
    model.result("pg3").feature("tblp1").set("probetag", "bnd1");
    model.result().export().create("tbl1", "Table");

    model.component("comp1").probe("bnd1").genResult(null);

    model.result("pg12").tag("pg3");

    model.study("std1").feature("param").set("pname", new String[]{"o0r1", "xoffs", "yoffs"});
    model.study("std1").feature("param")
         .set("plistarr", new String[]{redoxdirections, xpositions, ypositions});//<MOD/>
    model.study("std1").feature("param").set("punit", new String[]{"", "m", "m"});

    model.sol("sol1").attach("std1");
    model.sol("sol1").feature("st1").label("Compile Equations: Stationary");
    model.sol("sol1").feature("v1").label("Dependent Variables 1.1");
    model.sol("sol1").feature("v1").set("clistctrl", new String[]{"p1", "p1", "p1"});
    model.sol("sol1").feature("v1").set("cname", new String[]{"o0r1", "xoffs", "yoffs"});
    model.sol("sol1").feature("v1")
         .set("clist", new String[]{redoxdirections, xpositions, ypositions});//<MOD/>
    model.sol("sol1").feature("s1").label("Stationary Solver 1.1");
    model.sol("sol1").feature("s1").set("probesel", "none");
    model.sol("sol1").feature("s1").feature("dDef").label("Direct 2");
    model.sol("sol1").feature("s1").feature("aDef").label("Advanced 1");
    model.sol("sol1").feature("s1").feature("aDef").set("cachepattern", true);
    model.sol("sol1").feature("s1").feature("p1").label("Parametric 1.1");
    model.sol("sol1").feature("s1").feature("p1").set("control", "param");
    model.sol("sol1").feature("s1").feature("p1").set("pname", new String[]{"o0r1", "xoffs", "yoffs"});
    model.sol("sol1").feature("s1").feature("p1")
         .set("plistarr", new String[]{redoxdirections, xpositions, ypositions});//<MOD/>
    model.sol("sol1").feature("s1").feature("p1").set("punit", new String[]{"", "m", "m"});
    model.sol("sol1").feature("s1").feature("p1").set("pcontinuationmode", "no");
    model.sol("sol1").feature("s1").feature("fc1").label("Fully Coupled 1.1");
    model.sol("sol1").feature("s1").feature("fc1").set("linsolver", "i1");
    model.sol("sol1").feature("s1").feature("fc1").set("initstep", 0.01);
    model.sol("sol1").feature("s1").feature("fc1").set("minstep", 1.0E-6);
    model.sol("sol1").feature("s1").feature("fc1").set("maxiter", 50);
    model.sol("sol1").feature("s1").feature("i1").label("AMG, concentrations (tds)");
    model.sol("sol1").feature("s1").feature("i1").set("nlinnormuse", true);
    model.sol("sol1").feature("s1").feature("i1").set("maxlinit", 1000);
    model.sol("sol1").feature("s1").feature("i1").feature("ilDef").label("Incomplete LU 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").label("Multigrid 1.1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").set("prefun", "saamg");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").set("maxcoarsedof", 50000);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").set("saamgcompwise", true);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").set("usesmooth", false);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").label("Presmoother 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("soDef").label("SOR 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1").label("SOR Line 1.1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1")
         .set("linesweeptype", "ssor");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1").set("iter", 1);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1").set("linerelax", 0.7);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1").set("relax", 0.5);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").label("Postsmoother 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("soDef").label("SOR 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1").label("SOR Line 1.1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1")
         .set("linesweeptype", "ssor");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1").set("iter", 1);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1").set("linerelax", 0.7);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1").set("relax", 0.5);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").label("Coarse Solver 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").feature("dDef").label("Direct 2");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").feature("d1").label("Direct 1.1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").feature("d1")
         .set("linsolver", "pardiso");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").feature("d1")
         .set("pivotperturb", 1.0E-13);
    model.sol("sol1").feature("s1").feature("d1").label("Direct, concentrations (tds)");
    model.sol("sol1").feature("s1").feature("d1").set("linsolver", "pardiso");
    model.sol("sol1").feature("s1").feature("d1").set("pivotperturb", 1.0E-13);
    //model.sol("sol1").runAll();

    model.result().dataset("dset2").label("Probe Solution 2");
    model.result("pg6").set("view", "view2");
    model.result("pg6").feature("tbls1").set("dataformat", "filledtable");
    model.result("pg6").feature("surf1").set("resolution", "normal");
    model.result("pg7").set("view", "view3");
    model.result("pg7").set("xlabel", "X_Norm");
    model.result("pg7").set("ylabel", "Y_Norm");
    model.result("pg7").set("xlabelactive", false);
    model.result("pg7").set("ylabelactive", false);
    model.result("pg7").feature("tbls1").set("dataformat", "filledtable");
    model.result("pg8").label("Concentration, Cox, Streamline (tds)");
    model.result("pg8").set("titletype", "custom");
    model.result("pg8").set("prefixintitle", "Species Cox:");
    model.result("pg8").feature("str1").set("posmethod", "start");
    model.result("pg8").feature("str1").set("pointtype", "arrow");
    model.result("pg8").feature("str1").set("arrowlength", "logarithmic");
    model.result("pg8").feature("str1").set("resolution", "normal");
    model.result("pg8").feature("str1").feature("col").set("descr", "Concentration");
    model.result("pg8").feature("str1").feature("col").set("titletype", "custom");
    model.result("pg9").label("Concentration, Cox, Surface (tds)");
    model.result("pg9").set("titletype", "custom");
    model.result("pg9").set("prefixintitle", "Species Cox:");
    model.result("pg9").set("typeintitle", false);
    model.result("pg9").feature("surf1").set("descr", "Concentration");
    model.result("pg9").feature("surf1").set("resolution", "normal");
    model.result("pg10").label("Concentration, Cred, Streamline (tds)");
    model.result("pg10").set("titletype", "custom");
    model.result("pg10").set("prefixintitle", "Species Cred:");
    model.result("pg10").feature("str1")
         .set("expr", new String[]{"tds.tflux_Credx", "tds.tflux_Credy", "tds.tflux_Credz"});
    model.result("pg10").feature("str1").set("posmethod", "start");
    model.result("pg10").feature("str1").set("pointtype", "arrow");
    model.result("pg10").feature("str1").set("arrowlength", "logarithmic");
    model.result("pg10").feature("str1").set("resolution", "normal");
    model.result("pg10").feature("str1").feature("col").set("titletype", "custom");
    model.result("pg11").label("Concentration, Cred, Surface (tds)");
    model.result("pg11").set("titletype", "custom");
    model.result("pg11").set("prefixintitle", "Species Cred:");
    model.result("pg11").set("typeintitle", false);
    model.result("pg11").feature("surf1").set("resolution", "normal");
    model.result("pg3").label("Probe Plot Group 3");
    model.result("pg3").set("xlabel", "xoffs (m)");
    model.result("pg3").set("ylabel", "Current at UME (A), UME_FLUX");
    model.result("pg3").set("xlabelactive", false);
    model.result("pg3").set("ylabelactive", false);
    model.result().export("tbl1").set("filename", "data.txt");
    model.result().export("tbl1").set("ifexists", "append");

    model.component("comp1").probe("bnd1").genResult("none");

    model.sol("sol1").runAll();

    model.result().export("tbl1").run();

    return model;
  }
    
    public static Model runk(String reactivity_data_filename, double[] kreds, double[] koxs, double l, double xcoord, double ycoord, int o0r1){
        Model model = ModelUtil.create("Model");

    //<MOD>
	String krs = toString(kreds);
	String kos = toString(koxs);
	
	System.out.println("[DEBUG] x: " + xcoord + "[m] y: " + ycoord + "[m] krs: " + krs);
	System.out.println("[DEBUG] x: " + xcoord + "[m] y: " + ycoord + "[m] kos: " + kos);
	
	String cwd = getCWD();

    model.modelPath(cwd);
	//</MOD>

    model.label("SECM_BV.mph");

    model.param().set("a", "UME_Diam / 2", "UME RADIUS");
	//<MOD>
    model.param().set("Rg", RG + "*a", "Normalized Glass RADIUS");
    model.param().set("UME_Diam", UME_DIAM + "[um]", "UME DIAMETER");//</MOD
    model.param().set("X_UME", "Max_XY/2", "X-position of UME");
    model.param().set("Y_UME", "Max_XY/2", "Y-position of UME");
    //<MOD>
    model.param().set("Z_UME", l + "*a", "Z-position of UME");
    model.param().set("X_ADJ", "0[m]", "X-adjustment from the fit");
    model.param().set("Y_ADJ", "0[m]", "Y-adjustment from the fit");//</MOD
    model.param().set("Max_XY", "a*350", "Maximum x,y bound");
    model.param().set("Max_Z", "Max_XY *1.0", "Maximum z bound");
    model.param().set("L_UME", "Max_Z - Z_UME", "Length of UME");
	//<MOD>
    model.param().set("Cox_bulk", C_OX_BULK + "[mM]", "Bulk concentration of oxidized species");
    model.param().set("Cred_bulk", C_RED_BULK + "[mM]", "Bulk concentration of reduced species");
    model.param().set("Dox", D_OX + "[m^2/s]", "Diffusion coefficient - oxidised");
	model.param().set("Dred", D_RED + "[m^2/s]", "Diffusion coefficient - reduced");//</MOD
    model.param().set("k_UME", "1 [m/s]", "Rate constant at UME");
    model.param().set("k_PS", "10^logk [m/s]", "Rate constant at point source");
    model.param()
         .set("MaxMeshRadial_frac", "0.04", "Fraction of the radial features that will be the maximum mesh size");
    model.param().set("sl", "201", "grid side length");
    model.param().set("gridres", "0.1*a", "Grid resolution");
    model.param().set("logk", "-4.005", "base-10 logarithm of k/1m/s");
    //<MOD>
	model.param().set("xoffs", xcoord + "[m]");
    model.param().set("yoffs", ycoord + "[m]");//</MOD>
    model.param().set("factor", "10^(logk)/(1E-4)");
    model.param().set("kox_factor", "1");
    model.param().set("kred_factor", "1");

    model.component().create("comp1", true);

    model.component("comp1").geom().create("geom1", 3);

    model.component("comp1").curvedInterior(false);

    model.result().table().create("tbl1", "Table");

    model.component("comp1").func().create("int1", "Interpolation");
    model.component("comp1").func().create("an1", "Analytic");
    model.component("comp1").func().create("an2", "Analytic");
    model.component("comp1").func("int1").set("source", "file");
    model.component("comp1").func("int1").label("rate");
    model.component("comp1").func("int1").set("funcs", new String[][]{{"red_rate_map", "1"}, {"ox_rate_map", "2"}});
    //<MOD>
    model.component("comp1").func("int1").set("filename", cwd + "/" + reactivity_data_filename);//</MOD>
    model.component("comp1").func("int1").set("nargs", 2);
    model.component("comp1").func("int1").set("interp", "neighbor");
    model.component("comp1").func("int1").set("extrap", "value");
    model.component("comp1").func("int1").set("fununit", new String[]{"m/s", "m/s"});
    model.component("comp1").func("int1").set("argunit", new String[]{"m", "m"});
    model.component("comp1").func("an1").label("ShiftedKred");
    model.component("comp1").func("an1").set("funcname", "kred_func");
    model.component("comp1").func("an1")
         .set("expr", "red_rate_map(x+xoffs-0.5*Max_XY+X_ADJ,y+yoffs-0.5*Max_XY+Y_ADJ)*kred_factor");
    model.component("comp1").func("an1").set("args", new String[]{"x", "y"});
    model.component("comp1").func("an1").set("fununit", "m/s");
    model.component("comp1").func("an1").set("argunit", new String[]{"m", "m"});
    model.component("comp1").func("an1")
         .set("plotargs", new String[][]{{"x", "Max_XY*0.4", "Max_XY*0.6"}, {"y", "Max_XY*0.4", "Max_XY*0.6"}});
    model.component("comp1").func("an2").label("ShiftedKox");
    model.component("comp1").func("an2").set("funcname", "kox_func");
    model.component("comp1").func("an2")
         .set("expr", "ox_rate_map(x+xoffs-0.5*Max_XY+X_ADJ,y+yoffs-0.5*Max_XY+Y_ADJ)*kox_factor");
    model.component("comp1").func("an2").set("args", new String[]{"x", "y"});
    model.component("comp1").func("an2").set("fununit", "m/s");
    model.component("comp1").func("an2").set("argunit", new String[]{"m", "m"});
    model.component("comp1").func("an2")
         .set("plotargs", new String[][]{{"x", "Max_XY*0.4", "Max_XY*0.6"}, {"y", "Max_XY*0.4", "Max_XY*0.6"}});

    model.component("comp1").mesh().create("mesh1");

    model.component("comp1").geom("geom1").label("Simulation cell");
    model.component("comp1").geom("geom1").create("blk1", "Block");
    model.component("comp1").geom("geom1").feature("blk1").label("Bath");
    model.component("comp1").geom("geom1").feature("blk1").set("size", new String[]{"Max_XY", "Max_XY", "Max_Z"});
    model.component("comp1").geom("geom1").create("cyl1", "Cylinder");
    model.component("comp1").geom("geom1").feature("cyl1").label("UME");
    model.component("comp1").geom("geom1").feature("cyl1").set("pos", new String[]{"X_UME", "Y_UME", "Z_UME"});
    model.component("comp1").geom("geom1").feature("cyl1").set("r", "a");
    model.component("comp1").geom("geom1").feature("cyl1").set("h", "L_UME");
    model.component("comp1").geom("geom1").create("cyl3", "Cylinder");
    model.component("comp1").geom("geom1").feature("cyl3").label("UME_Ins");
    model.component("comp1").geom("geom1").feature("cyl3").set("pos", new String[]{"X_UME", "Y_UME", "Z_UME"});
    model.component("comp1").geom("geom1").feature("cyl3").set("r", "Rg");
    model.component("comp1").geom("geom1").feature("cyl3").set("h", "L_UME");
    model.component("comp1").geom("geom1").create("blk2", "Block");
    model.component("comp1").geom("geom1").feature("blk2").active(false);
    model.component("comp1").geom("geom1").feature("blk2").label("box00");
    model.component("comp1").geom("geom1").feature("blk2").set("pos", new String[]{"42.25*a", "43.75*a", "-1E-5"});
    model.component("comp1").geom("geom1").feature("blk2").set("size", new String[]{"0.5*a", "0.5*a", "1E-5"});
    model.component("comp1").geom("geom1").run();
    model.component("comp1").geom("geom1").run("fin");

    model.view().create("view2", 2);
    model.view().create("view3", 2);

    model.component("comp1").physics().create("tds", "DilutedSpecies", "geom1");
    model.component("comp1").physics("tds").field("concentration").field("Cox");
    model.component("comp1").physics("tds").field("concentration").component(new String[]{"Cox", "Cred"});
    model.component("comp1").physics("tds").selection().set(1);
    model.component("comp1").physics("tds").create("fl1", "FluxBoundary", 2);
    model.component("comp1").physics("tds").feature("fl1").selection().set(12);
    model.component("comp1").physics("tds").create("fl2", "FluxBoundary", 2);
    model.component("comp1").physics("tds").feature("fl2").selection().set(3);
    model.component("comp1").physics("tds").create("conc1", "Concentration", 2);
    model.component("comp1").physics("tds").feature("conc1").selection().set(1, 2, 4, 5, 18);

    model.component("comp1").mesh("mesh1").create("ftet1", "FreeTet");
    model.component("comp1").mesh("mesh1").feature("ftet1").selection().geom("geom1", 3);
    model.component("comp1").mesh("mesh1").feature("ftet1").selection().set(1);
    model.component("comp1").mesh("mesh1").feature("ftet1").create("size1", "Size");
    model.component("comp1").mesh("mesh1").feature("ftet1").create("size2", "Size");
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size1").selection().geom("geom1", 2);
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size1").selection().set(12);
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size2").selection().geom("geom1", 2);
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size2").selection().set(3);

    model.component("comp1").probe().create("bnd1", "Boundary");
    model.component("comp1").probe("bnd1").selection().set(12);

    model.result().table("tbl1").label("Probe Table 1");

    model.component("comp1").view("view1").set("transparency", true);
    model.view("view2").axis().set("xmin", 397.5);
    model.view("view2").axis().set("xmax", 865);
    model.view("view2").axis().set("ymin", 494.375);
    model.view("view2").axis().set("ymax", 755.625);
    model.view("view2").axis().set("viewscaletype", "automatic");
    model.view("view3").axis().set("xmin", 41.69999694824219);
    model.view("view3").axis().set("xmax", 59.30000305175781);
    model.view("view3").axis().set("ymin", 42.8861083984375);
    model.view("view3").axis().set("ymax", 57.1138916015625);
    model.view("view3").axis().set("viewscaletype", "automatic");

    model.component("comp1").physics("tds").prop("EquationForm").set("form", "Stationary");
    model.component("comp1").physics("tds").prop("TransportMechanism").set("Convection", false);
    model.component("comp1").physics("tds").feature("cdm1")
         .set("D_Cox", new String[][]{{"Dox"}, {"0"}, {"0"}, {"0"}, {"Dox"}, {"0"}, {"0"}, {"0"}, {"Dox"}});
    model.component("comp1").physics("tds").feature("cdm1")
         .set("D_Cred", new String[][]{{"Dred"}, {"0"}, {"0"}, {"0"}, {"Dred"}, {"0"}, {"0"}, {"0"}, {"Dred"}});
    model.component("comp1").physics("tds").feature("init1")
         .set("initc", new String[][]{{"Cox_bulk"}, {"Cred_bulk"}});
    model.component("comp1").physics("tds").feature("fl1").set("species", new int[][]{{1}, {1}});
    if(o0r1 == 0){
		model.component("comp1").physics("tds").feature("fl1")
         .set("J0", new String[][]{{"+k_UME*Cred"}, {"-k_UME*Cred"}});
	}
	else{
		model.component("comp1").physics("tds").feature("fl1")
         .set("J0", new String[][]{{"-k_UME*Cox"}, {"+k_UME*Cox"}});
	}
    model.component("comp1").physics("tds").feature("fl1").label("Flux_UME");
    model.component("comp1").physics("tds").feature("fl2").set("species", new int[][]{{1}, {1}});
    model.component("comp1").physics("tds").feature("fl2")
         .set("J0", new String[][]{{"(kox_func(x,y) * Cred - kred_func(x,y) * Cox)"}, {"(kred_func(x,y) * Cox - kox_func(x,y) * Cred)"}});
    model.component("comp1").physics("tds").feature("fl2").label("Flux_PS");
    model.component("comp1").physics("tds").feature("conc1").set("species", new int[][]{{1}, {1}});
    model.component("comp1").physics("tds").feature("conc1").set("c0", new String[][]{{"Cox_bulk"}, {"Cred_bulk"}});
    model.component("comp1").physics("tds").feature("conc1").label("Bulk");

    model.component("comp1").mesh("mesh1").feature("size").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("size").set("hmax", "80[um]");
    model.component("comp1").mesh("mesh1").feature("size").set("hmin", "10[um]");
    model.component("comp1").mesh("mesh1").feature("size").set("hgrad", 1.15);
    model.component("comp1").mesh("mesh1").feature("ftet1").set("optcurved", false);
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size1").label("UME");
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size1").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size1")
         .set("hmax", "a*MaxMeshRadial_frac*0.18");
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size1").set("hmaxactive", true);
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size1").set("hmin", 1.8E-5);
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size1").set("hminactive", false);
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size2").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size2").set("hgrad", "1.10");
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size2").set("hgradactive", true);
    model.component("comp1").mesh("mesh1").run();

    model.component("comp1").probe("bnd1").label("UME_FLUX");
    model.component("comp1").probe("bnd1").set("type", "integral");
    model.component("comp1").probe("bnd1").set("probename", "I_UME");
    model.component("comp1").probe("bnd1").set("expr", "tds.ntflux_Cred*F_const");
    model.component("comp1").probe("bnd1").set("unit", "A");
    model.component("comp1").probe("bnd1").set("descractive", true);
    model.component("comp1").probe("bnd1").set("descr", "Current at UME");
    model.component("comp1").probe("bnd1").set("table", "tbl1");
    model.component("comp1").probe("bnd1").set("window", "window1");

    model.study().create("std1");
    model.study("std1").create("param", "Parametric");
    model.study("std1").create("stat", "Stationary");

    model.sol().create("sol1");
    model.sol("sol1").study("std1");
    model.sol("sol1").attach("std1");
    model.sol("sol1").create("st1", "StudyStep");
    model.sol("sol1").create("v1", "Variables");
    model.sol("sol1").create("s1", "Stationary");
    model.sol("sol1").feature("s1").create("p1", "Parametric");
    model.sol("sol1").feature("s1").create("fc1", "FullyCoupled");
    model.sol("sol1").feature("s1").create("i1", "Iterative");
    model.sol("sol1").feature("s1").create("d1", "Direct");
    model.sol("sol1").feature("s1").feature("i1").create("mg1", "Multigrid");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").create("sl1", "SORLine");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").create("sl1", "SORLine");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").create("d1", "Direct");
    model.sol("sol1").feature("s1").feature().remove("fcDef");

    model.result().dataset().create("dset2", "Solution");
    model.result().dataset().create("int1", "Integral");
    model.result().dataset("dset2").set("probetag", "bnd1");
    model.result().dataset("int1").set("probetag", "bnd1");
    model.result().dataset("int1").set("data", "dset2");
    model.result().dataset("int1").selection().geom("geom1", 2);
    model.result().dataset("int1").selection().set(12);
    model.result().numerical().create("pev1", "EvalPoint");
    model.result().numerical("pev1").set("probetag", "bnd1");
    model.result().create("pg6", "PlotGroup2D");
    model.result().create("pg7", "PlotGroup2D");
    model.result().create("pg3", "PlotGroup1D");
    model.result().create("pg8", "PlotGroup3D");
    model.result().create("pg9", "PlotGroup3D");
    model.result().create("pg10", "PlotGroup3D");
    model.result().create("pg11", "PlotGroup3D");
    model.result("pg6").create("tbls1", "TableSurface");
    model.result("pg6").create("surf1", "Surface");
    model.result("pg6").feature("surf1").set("expr", "comp1.Cox");
    model.result("pg7").create("tbls1", "TableSurface");
    model.result("pg3").set("probetag", "window1");
    model.result("pg3").create("tblp1", "Table");
    model.result("pg3").feature("tblp1").set("probetag", "bnd1");
    model.result("pg8").create("str1", "Streamline");
    model.result("pg8").feature("str1").create("col", "Color");
    model.result("pg9").create("surf1", "Surface");
    model.result("pg10").create("str1", "Streamline");
    model.result("pg10").feature("str1").create("col", "Color");
    model.result("pg10").feature("str1").feature("col").set("expr", "Cred");
    model.result("pg11").create("surf1", "Surface");
    model.result("pg11").feature("surf1").set("expr", "Cred");
    model.result().export().create("tbl1", "Table");

    model.component("comp1").probe("bnd1").genResult(null);

    model.result("pg12").tag("pg3");

    model.study("std1").feature("param").set("pname", new String[]{"kox_factor", "kred_factor"});
	//<MOD>
    model.study("std1").feature("param").set("plistarr", new String[]{kos, krs});//</MOD>
    model.study("std1").feature("param").set("punit", new String[]{"", ""});

    model.sol("sol1").attach("std1");
    model.sol("sol1").feature("st1").label("Compile Equations: Stationary");
    model.sol("sol1").feature("v1").label("Dependent Variables 1.1");
    model.sol("sol1").feature("v1").set("clistctrl", new String[]{"p1", "p1"});
    model.sol("sol1").feature("v1").set("cname", new String[]{"kox_factor", "kred_factor"});
    //<MOD>
	model.sol("sol1").feature("v1").set("clist", new String[]{kos, krs});//</MOD>
    model.sol("sol1").feature("s1").label("Stationary Solver 1.1");
    model.sol("sol1").feature("s1").set("probesel", "none");
    model.sol("sol1").feature("s1").feature("dDef").label("Direct 2");
    model.sol("sol1").feature("s1").feature("aDef").label("Advanced 1");
    model.sol("sol1").feature("s1").feature("aDef").set("cachepattern", true);
    model.sol("sol1").feature("s1").feature("p1").label("Parametric 1.1");
    model.sol("sol1").feature("s1").feature("p1").set("control", "param");
    model.sol("sol1").feature("s1").feature("p1").set("pname", new String[]{"kox_factor", "kred_factor"});
	//<MOD>
    model.sol("sol1").feature("s1").feature("p1").set("plistarr", new String[]{kos, krs});//</MOD>
    model.sol("sol1").feature("s1").feature("p1").set("punit", new String[]{"", ""});
    model.sol("sol1").feature("s1").feature("p1").set("pcontinuationmode", "no");
    model.sol("sol1").feature("s1").feature("fc1").label("Fully Coupled 1.1");
    model.sol("sol1").feature("s1").feature("fc1").set("linsolver", "i1");
    model.sol("sol1").feature("s1").feature("fc1").set("initstep", 0.01);
    model.sol("sol1").feature("s1").feature("fc1").set("minstep", 1.0E-6);
    model.sol("sol1").feature("s1").feature("fc1").set("maxiter", 50);
    model.sol("sol1").feature("s1").feature("i1").label("AMG, concentrations (tds)");
    model.sol("sol1").feature("s1").feature("i1").set("nlinnormuse", true);
    model.sol("sol1").feature("s1").feature("i1").set("maxlinit", 1000);
    model.sol("sol1").feature("s1").feature("i1").feature("ilDef").label("Incomplete LU 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").label("Multigrid 1.1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").set("prefun", "saamg");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").set("maxcoarsedof", 50000);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").set("saamgcompwise", true);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").set("usesmooth", false);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").label("Presmoother 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("soDef").label("SOR 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1").label("SOR Line 1.1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1")
         .set("linesweeptype", "ssor");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1").set("iter", 1);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1").set("linerelax", 0.7);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1").set("relax", 0.5);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").label("Postsmoother 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("soDef").label("SOR 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1").label("SOR Line 1.1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1")
         .set("linesweeptype", "ssor");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1").set("iter", 1);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1").set("linerelax", 0.7);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1").set("relax", 0.5);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").label("Coarse Solver 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").feature("dDef").label("Direct 2");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").feature("d1").label("Direct 1.1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").feature("d1")
         .set("linsolver", "pardiso");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").feature("d1")
         .set("pivotperturb", 1.0E-13);
    model.sol("sol1").feature("s1").feature("d1").label("Direct, concentrations (tds)");
    model.sol("sol1").feature("s1").feature("d1").set("linsolver", "pardiso");
    model.sol("sol1").feature("s1").feature("d1").set("pivotperturb", 1.0E-13);
    model.sol("sol1").runAll();

    model.result().dataset("dset2").label("Probe Solution 2");
    model.result("pg6").set("view", "view2");
    model.result("pg6").feature("tbls1").set("dataformat", "filledtable");
    model.result("pg6").feature("surf1").set("resolution", "normal");
    model.result("pg7").set("view", "view3");
    model.result("pg7").set("xlabel", "X_Norm");
    model.result("pg7").set("ylabel", "Y_Norm");
    model.result("pg7").set("xlabelactive", false);
    model.result("pg7").set("ylabelactive", false);
    model.result("pg7").feature("tbls1").set("dataformat", "filledtable");
    model.result("pg3").label("Probe Plot Group 3");
    model.result("pg3").set("xlabel", "xoffs (m)");
    model.result("pg3").set("ylabel", "Current at UME (A), UME_FLUX");
    model.result("pg3").set("windowtitle", "Probe Plot 2");
    model.result("pg3").set("xlabelactive", false);
    model.result("pg3").set("ylabelactive", false);
    model.result("pg8").label("Concentration, Cox, Streamline (tds)");
    model.result("pg8").set("titletype", "custom");
    model.result("pg8").set("prefixintitle", "Species Cox:");
    model.result("pg8").feature("str1").set("posmethod", "start");
    model.result("pg8").feature("str1").set("pointtype", "arrow");
    model.result("pg8").feature("str1").set("arrowlength", "logarithmic");
    model.result("pg8").feature("str1").set("resolution", "normal");
    model.result("pg8").feature("str1").feature("col").set("descr", "Concentration");
    model.result("pg8").feature("str1").feature("col").set("titletype", "custom");
    model.result("pg9").label("Concentration, Cox, Surface (tds)");
    model.result("pg9").set("titletype", "custom");
    model.result("pg9").set("prefixintitle", "Species Cox:");
    model.result("pg9").set("typeintitle", false);
    model.result("pg9").feature("surf1").set("descr", "Concentration");
    model.result("pg9").feature("surf1").set("resolution", "normal");
    model.result("pg10").label("Concentration, Cred, Streamline (tds)");
    model.result("pg10").set("titletype", "custom");
    model.result("pg10").set("prefixintitle", "Species Cred:");
    model.result("pg10").feature("str1")
         .set("expr", new String[]{"tds.tflux_Credx", "tds.tflux_Credy", "tds.tflux_Credz"});
    model.result("pg10").feature("str1").set("posmethod", "start");
    model.result("pg10").feature("str1").set("pointtype", "arrow");
    model.result("pg10").feature("str1").set("arrowlength", "logarithmic");
    model.result("pg10").feature("str1").set("resolution", "normal");
    model.result("pg10").feature("str1").feature("col").set("titletype", "custom");
    model.result("pg11").label("Concentration, Cred, Surface (tds)");
    model.result("pg11").set("titletype", "custom");
    model.result("pg11").set("prefixintitle", "Species Cred:");
    model.result("pg11").set("typeintitle", false);
    model.result("pg11").feature("surf1").set("resolution", "normal");
    model.result().export("tbl1").set("filename", "data.txt");
    model.result().export("tbl1").set("ifexists", "append");

    model.component("comp1").probe("bnd1").genResult("none");

    model.sol("sol1").runAll();

    model.result().export("tbl1").run();

    return model;
    }
    
    public static void main(String[] args) throws NumberFormatException, FileNotFoundException {
        
        
        try{
            String control_file = "if_252mV-3.csv";
            System.out.println(getDateStamp() + " Reading in control file: " + control_file);
            readSECMInfo(control_file);
			
            System.out.println(getDateStamp() + " Finding initial logk");
            double initial_r = 0;
            double initial_l = 1.0;			
			double initial_e = 0.0;
			double initial_logk = findFirstLogK(initial_l, initial_e, true);
			System.out.println("initialk: " + initial_logk);
			
			list_l = new LinkedList<Double>();
			list_e = new LinkedList<Double>();
			list_logk = new LinkedList<Double>();
			list_gridsum = new LinkedList<Integer>();
			list_data = new LinkedList<Double[]>();
			list_r = new LinkedList<Double>();
			list_xa = new LinkedList<Double>();
			list_ya = new LinkedList<Double>();

			int result = runFit(control_file, initial_l, initial_e, initial_logk, initial_r, true);
		  	logEndCondition(result);
			
        }
        catch(Exception e){
            System.out.println(e.toString());
			System.out.println(e.getMessage());
            e.printStackTrace();
        }
    }
    
    /*
    Fitting methods and fields
    */
    /**
     * Runs the fitting procedure using Levenberg-Marquardt
     * @param filename The name of the file containing the data to be fit to (for logging purposes). {@link #readSECMInfo(java.lang.String)} should be used before this method to load-in the necessary data.
     * @param firstl The initial guess for L
     * @param firstlogk The initial guess for logk. Use of {@link #findFirstLogK(double, boolean)} to obtain this value is recommended.
     * @param verbose flag to be ticked if extra data is desired (the current and derivatives at each iteration)
     * @return The status of the fitting procedure. 
     * <p>{@link #EXECUTED_OK} if the process converged,</p>
     * <p>{@link #MAX_ITERATIONS_REACHED} if there is no convergence after {@link #MAX_ITERATIONS} iterations,</p>
     * <p>or {@link #MAX_LAMBDA_REACHED} if lambda exceeds {@link #MAX_LAMBDA} when trying to compute the next iteration's parameters.</p>
     * @throws FileNotFoundException
     * @throws NumberFormatException
     * @throws IOException 
     */
    static int runFit(String filename, double firstl, double firste, double firstlogk, double firstr, boolean verbose) throws FileNotFoundException, NumberFormatException, IOException{
        String[] param_labels = new String[]{"L", "e", "logk", "dilation/erosion", "x-adjustment", "y-adjustment"};
        double[] experimental = true_image;
        double lambda = 0;
        double ssr = 0;
        double last_ssr;
        double first_xa = 0;
        double first_ya = 0;
        double l = firstl;
        double e = firste;
        double logk = firstlogk;
        double r = firstr;
        double xa = first_xa;
        double ya = first_ya;
        double last_l;
        double last_e;
        double last_logk;
        double last_r;
        double last_xa = xa;
        double last_ya = ya;
        
        //first iteration
        //compute the currents and residuals for the initial parameter guesses
        double[] curr = runModel(firstl, firste, firstlogk, firstr, first_xa, first_ya);
        double[] residuals = subtract(experimental, curr);
        
        //compute the derivatives of the current with respect to the parameters
        double[] curr_dl = runModel(firstl + L_PERTURB, firste, firstlogk, firstr, first_xa, first_ya);
        double[] dl = multiply(subtract(curr_dl, curr), 1.0/L_PERTURB);
        
        double[] curr_de = runModel(firstl, firste + E_PERTURB, firstlogk, firstr, first_xa, first_ya);
        double[] de = multiply(subtract(curr_de, curr), 1.0/E_PERTURB);
        
        double[] curr_dlogk = runModel(firstl, firste, firstlogk + LOGK_PERTURB, firstr, first_xa, first_ya);
        double[] dlogk = multiply(subtract(curr_dlogk, curr), 1.0/LOGK_PERTURB);
        
        double[] curr_dr = runModel(firstl, firste, firstlogk, firstr + R_PERTURB, first_xa, first_ya);
        double[] dr = multiply(subtract(curr_dr, curr), 1.0/R_PERTURB);
        
        double[] curr_dxa = runModel(firstl, firste, firstlogk, firstr, first_xa + ADJ_X_PERTURB, first_ya);
        double[] dxa = multiply(subtract(curr_dxa, curr), 1.0/ADJ_X_PERTURB);
        
        double[] curr_dya = runModel(firstl, firste, firstlogk, firstr, first_xa, first_ya + ADJ_Y_PERTURB);
        double[] dya = multiply(subtract(curr_dya, curr), 1.0/ADJ_Y_PERTURB);
        
        //construct the Jacobian, and perform an iteration of Levenberg-Marquardt
        //[dl dk dr]
        double[][] J = appendColumn(appendColumn(appendColumn(appendColumn(appendColumn(dl, de), dlogk), dr), dxa), dya);
        double[][] JT = transpose(J);
        double[][] JTJ = multiply(JT, J);
        double[][] DTD = new double[JTJ.length][JTJ.length];
        for(int i = 0; i < DTD.length; i++){
            for(int j = 0; j < DTD.length; j++){
                if(i != j){
                    DTD[i][j] = 0;
                }
                else{
                    DTD[i][j] = JTJ[i][j];
                }
            }
        }
        double[][] lam_DTD = multiply(DTD, lambda);
        double[][] JTJinv = invert(add(JTJ, lam_DTD));
        double[] delta_c = multiply(multiply(JTJinv, JT), residuals);
        ssr = sumSquare(residuals);
        //initialize 'previous iteration' data
        last_ssr = ssr;
        last_l = l;
        last_e = e;
        last_logk = logk;
        last_r = r;
		last_xa = xa;
		last_ya = ya;
		logInitialGuesses(filename, param_labels, new double[]{l, e, logk, r, xa, ya}, ssr);
        l = applyLLimit(last_l + round(delta_c[0], L_DECIMALS));
        e = applyELimit(last_e + round(delta_c[1], E_DECIMALS));
        logk = applyLOGKLimit(last_logk + round(delta_c[2], LOGK_DECIMALS));
        r = applyRLimit(last_r + round(delta_c[3], R_DECIMALS));
        xa = applyADJLimit(last_xa + round(delta_c[4], ADJ_DECIMALS));
        ya = applyADJLimit(last_ya + round(delta_c[5], ADJ_DECIMALS));
        
        boolean converged = false;
        int iterations = 1;
        if(verbose){
            writeIteration("Iteration_1.txt", physical_xs, physical_ys, me_o0r1, curr, dl, de, dlogk, dr, dxa, dya, residuals);
        }
        //subsequent iterations
        while(!converged && iterations <= MAX_ITERATIONS){
            iterations ++;
            System.out.println("Iteration " + iterations); 
            System.out.println("Lambda 0");
            //First lambda
            curr = runModel(l, e, logk, r, xa, ya);
            residuals = subtract(experimental, curr);
            ssr = sumSquare(residuals);
            boolean lambda_ok = ssr < last_ssr;
            logIteration(iterations, new double[]{DTD[0][0], DTD[1][1], DTD[2][2], DTD[3][3], DTD[4][4], DTD[5][5]}, lambda, param_labels, new double[]{l, e, logk, r, xa, ya}, ssr, lambda_ok);

            //subsequent lambdas (if necessary)
            while(!lambda_ok && lambda <= MAX_LAMBDA){
                lambda = nextLambda(lambda);
                System.out.println("Lambda " + lambda);
                lam_DTD = multiply(DTD, lambda);
                JTJinv = invert(add(JTJ, lam_DTD));
                delta_c = multiply(multiply(JTJinv, JT), residuals);
                l = applyLLimit(last_l + round(delta_c[0], L_DECIMALS));
                e = applyELimit(last_e + round(delta_c[1], E_DECIMALS));
                logk = applyLOGKLimit(last_logk + round(delta_c[2], LOGK_DECIMALS));
                r = applyRLimit(last_r + round(delta_c[3], R_DECIMALS));
                xa = applyADJLimit(last_xa + round(delta_c[4], ADJ_DECIMALS));
                ya = applyADJLimit(last_ya + round(delta_c[5], ADJ_DECIMALS));
                curr = runModel(l, e, logk, r, xa, ya);
                residuals = subtract(experimental, curr);
                ssr = sumSquare(residuals);
                lambda_ok = ssr < last_ssr;
                logLambda(lambda, param_labels, new double[]{l, e, logk, r, xa, ya}, ssr, lambda_ok);
                //if the prescribed changes to the parameters are small enough, declare convergence
                if(Math.abs(delta_c[0]) < 0.5*Math.pow(10, -L_DECIMALS) && Math.abs(delta_c[1]) < 0.5*Math.pow(10, -E_DECIMALS) && Math.abs(delta_c[2]) < 0.5*Math.pow(10, -LOGK_DECIMALS) && Math.abs(delta_c[3]) < 0.5*Math.pow(10, -R_DECIMALS) && Math.abs(delta_c[4]) < 0.5*Math.pow(10, -ADJ_DECIMALS) && Math.abs(delta_c[5]) < 0.5*Math.pow(10, -ADJ_DECIMALS)){
                    System.out.println("DeltaC: " + delta_c[0] + "\t" + delta_c[1] + "\t" + delta_c[2] + "\t" + delta_c[3] + "\t" + delta_c[4] + "\t" + delta_c[5]);
                    lambda_ok = true;
                    converged = true;
                    return EXECUTED_OK;
                }
            }
            //check if the loop ended due to lambda hitting its maximum
            if(!lambda_ok){
                return MAX_LAMBDA_REACHED;
            }
            //end the process if this was the last iteration.
            if(iterations >= MAX_ITERATIONS){
                return MAX_ITERATIONS_REACHED;
            }
            //Experimental: (use the set of parameters with the lowest residuals to-date)
            int lowest_sim = findLowestSSR();
                if(lowest_sim < list_l.size() - 1){
                    l = list_l.get(lowest_sim);
                    e = list_e.get(lowest_sim);
                    logk = list_logk.get(lowest_sim);
                    r = list_r.get(lowest_sim);
                    xa = list_xa.get(lowest_sim);
                    ya = list_ya.get(lowest_sim);
                    curr = convertToRegularDouble(list_data.get(lowest_sim));
                    residuals = subtract(experimental, curr);
                    ssr = sumSquare(residuals);
                    lambda_ok = true;
                    System.out.println(getDateStamp() + ": found lower: L: " + l + "; e: " + e + "; logk: " + logk + "; erosion/dilation: " + r + "; x-adjustment: " + xa + "; y-adjustment: " + ya);
                    logLambda(-1, param_labels, new double[]{l, e, logk, r, xa, ya}, ssr, lambda_ok);
            }
            //end of Experimental bit
            //update the previous iteration parameters
            last_ssr = ssr;
            last_l = l;
            last_e = e;
            last_logk = logk;
            last_r = r;
            last_xa = xa;
            last_ya = ya;
            int[][] de_grid = applyDilationErosion(r);
            int gridcount = getSum(de_grid);
            //compute or look up L derivative
            int query = checkList(l - L_PERTURB, e, logk, gridcount, xa, ya);
            if(query == -1){
                curr_dl = runModel(l + L_PERTURB, e, logk, r, xa, ya);
                dl = multiply(subtract(curr_dl, curr), 1.0/L_PERTURB);
            }
            else{
                curr_dl = runModel(l - L_PERTURB, e, logk, r, xa, ya);
                //NOTE: this curr_dl is effectively simulated as being perturbed negatively
                dl = multiply(subtract(curr_dl, curr), -1.0/L_PERTURB);
            }
            
            //compute or look up e derivative
            query = checkList(l, e - E_PERTURB, logk, gridcount, xa, ya);
            if(query == -1){
                curr_de = runModel(l, e + E_PERTURB, logk, r, xa, ya);
                de = multiply(subtract(curr_de, curr), 1.0/E_PERTURB);
            }
            else{
                curr_de = runModel(l, e - E_PERTURB, logk, r, xa, ya);
                //NOTE: this curr_de is effectively simulated as being perturbed negatively
                de = multiply(subtract(curr_de, curr), -1.0/E_PERTURB);
            }
            
            //compute or look up logk derivative
            query = checkList(l, e, logk - LOGK_PERTURB, gridcount, xa, ya);
            if(query == -1){
                curr_dlogk = runModel(l, e, logk + LOGK_PERTURB, r, xa, ya);
                dlogk = multiply(subtract(curr_dlogk, curr), 1.0/E_PERTURB);
            }
            else{
                curr_dlogk = runModel(l, e, logk - LOGK_PERTURB, r, xa, ya);
                //NOTE: this curr_dlogkk is effectively simulated as being perturbed negatively
                dlogk = multiply(subtract(curr_dlogk, curr), -1.0/LOGK_PERTURB);
            }
            
            //compute or look up xa derivative
            query = checkList(l, e, logk, gridcount, xa - ADJ_X_PERTURB, ya);
            if(query == -1){
                curr_dxa = runModel(l, e, logk, r, xa + ADJ_X_PERTURB, ya);
                dxa = multiply(subtract(curr_dxa, curr), 1.0/ADJ_X_PERTURB);
            }
            else{
                curr_dxa = runModel(l, e, logk, r, xa - ADJ_X_PERTURB, ya);
                dxa = multiply(subtract(curr_dxa, curr), -1.0/ADJ_X_PERTURB);
            }
            
            //compute or look up ya derivative
            query = checkList(l, e, logk, gridcount, xa, ya - ADJ_Y_PERTURB);
            if(query == -1){
                curr_dya = runModel(l, e, logk, r, xa, ya + ADJ_Y_PERTURB);
                dya = multiply(subtract(curr_dya, curr), 1.0/ADJ_Y_PERTURB);
            }
            else{
                curr_dya = runModel(l, e, logk, r, xa, ya - ADJ_Y_PERTURB);
                dya = multiply(subtract(curr_dya, curr), -1.0/ADJ_Y_PERTURB);
            }
            
            //compute or look up r derivative
            de_grid = applyDilationErosion(r - R_PERTURB);
            gridcount = getSum(de_grid);
            query = checkList(l, e, logk, gridcount, xa, ya);
            if(query == -1){
                curr_dr = runModel(l, e, logk, r + R_PERTURB, xa, ya);
                dr = multiply(subtract(curr_dr, curr), 1.0/R_PERTURB);
            }
            else{
                curr_dr = runModel(l, e, logk, r - R_PERTURB, xa, ya);
                //NOTE: this curr_dr is effectively simulated as being perturbed negatively
                dr = multiply(subtract(curr_dr, curr), -1.0/R_PERTURB);
            }
            

            if(verbose){
                writeIteration("Iteration_" + iterations + ".txt", physical_xs, physical_ys, me_o0r1, curr, dl, de, dlogk, dr, dxa, dya, residuals);
            }
            //construct the jacobian
            J = appendColumn(appendColumn(appendColumn(appendColumn(appendColumn(dl, de), dlogk), dr), dxa), dya);
            JT = transpose(J);
            JTJ = multiply(JT, J);
            for(int i = 0; i < DTD.length; i++){
                if(JTJ[i][i] > DTD[i][i]){
                    DTD[i][i] = JTJ[i][i];
                }
            }
            lambda = 0;
            lam_DTD = multiply(DTD, lambda);
                if(determinant(add(JTJ, lam_DTD)) != 0.0){
                    JTJinv = invert(add(JTJ, lam_DTD));
                    delta_c = multiply(multiply(JTJinv, JT), residuals);
                    l = applyLLimit(last_l + round(delta_c[0], L_DECIMALS));
                    e = applyELimit(last_e + round(delta_c[1], E_DECIMALS));
                    logk = applyLOGKLimit(last_logk + round(delta_c[2], LOGK_DECIMALS));
                    r = applyRLimit(last_r + round(delta_c[3], R_DECIMALS));
                    xa = applyADJLimit(last_xa + round(delta_c[4], ADJ_DECIMALS));
                    ya = applyADJLimit(last_ya + round(delta_c[5], ADJ_DECIMALS));
                    //if the prescribed changes to the parameters are small enough, declare convergence
                    if(Math.abs(delta_c[0]) < 0.5*Math.pow(10, -L_DECIMALS) && Math.abs(delta_c[1]) < 0.5*Math.pow(10, -E_DECIMALS) && Math.abs(delta_c[2]) < 0.5*Math.pow(10, -LOGK_DECIMALS) && Math.abs(delta_c[3]) < 0.5*Math.pow(10, -R_DECIMALS) && Math.abs(delta_c[4]) < 0.5*Math.pow(10, -ADJ_DECIMALS) && Math.abs(delta_c[5]) < 0.5*Math.pow(10, -ADJ_DECIMALS)){
                        System.out.println("DeltaC: " + delta_c[0] + "\t" + delta_c[1] + "\t" + delta_c[2] + "\t" + delta_c[3] + "\t" + delta_c[4] + "\t" + delta_c[5]);
                        converged = true;
                        return EXECUTED_OK;
                    }
                }
        }
        
        return MAX_ITERATIONS_REACHED;
    }
    
    /**
     * Checks if the parameters have already been used in a previous iteration, then either returns the previously calculated currents or simulates the currents.
     * @param l the L parameter
     * @param logk the base-ten logarithm of the k parameter
     * @return The currents of the microelectrode at the specified relative positions with respect to the reactive feature. 
     * The positions follow the same order that their currents are defined in the control file.
     * @throws FileNotFoundException 
     */
    static double[] runModel(double l, double e, double logk, double r, double xa, double ya) throws FileNotFoundException, IOException{
        int[][] de_grid = applyDilationErosion(r);
        int gridcount = getSum(de_grid);
        int index = checkList(l, e, logk, gridcount, xa, ya);
        if(index == -1){
            System.out.println(getDateStamp() + ":  simulating: L: " + l + "; e: " + e + "; logk: " + logk + "; erosion/dilation: " + r + "; x-adjustment: " + xa + "; y-adjustment: " + ya);
            double kred = Math.exp(-ALPHA*F*e)*Math.pow(10,logk);
            double kox = Math.exp((1.0-ALPHA)*F*e)*Math.pow(10,logk);
			
			if(kred > MAX_KRED){
				double factor = MAX_KRED / kred;
				kred *= factor;
				kox *= factor;
			}
			
			if(kox > MAX_KOX){
				double factor = MAX_KOX / kox;
				kred *= factor;
				kox *= factor;
			}
			
            writeReactivityFile(de_grid, kred, kox);
            Model temp = run(reactivity_mapfile, l, xa, ya, physical_xs, physical_ys, me_o0r1);
            double[] data = readData();
            eraseDataFile();
            addToList(l, e, logk, r, gridcount, xa, ya, data);
            return data;
        }
        else{
            System.out.println(getDateStamp() + ":      lookup: L: " + l + "; e: " + e + "; logk: " + logk + "; erosion/dilation: " + r + "; x-adjustment: " + xa + "; y-adjustment: " + ya);
            return convertToRegularDouble(list_data.get(index));
        }
    }
    
    /**
     * Prints a progress bar to the console
     * @param prog the progress so far
     * @param total the maximum possible value for progress
     * @return The progress bar as a string. Note: the first character of this string is a '\r'.
     */
    static String pBar(int prog, int total){
	int len = 25;
	int filled = (len*prog) / total;
	int percent = (100*prog) / total;
	String bar = " ";
	for(int i = 1; i <= filled; i++){
	  bar = bar.concat("\u2588");
	}
	for(int i = filled + 1; i <= len; i++){
	  bar = bar.concat("\u2591");
	}
	bar = bar.concat("  " + percent + "% complete.");
	return "\r " + bar;
    }
    
    /**
     * Returns The current ISO8601 timestamp (minus the timezone information)
     * @return The current time and date as yyyy-mm-ddThh:mm:ss
     */
    static String getDateStamp(){
        Calendar cl = Calendar.getInstance();
        int year = cl.get(Calendar.YEAR);
        int month = cl.get(Calendar.MONTH) + 1;
        int day = cl.get(Calendar.DAY_OF_MONTH);
        int hour = cl.get(Calendar.HOUR_OF_DAY);
        int minute = cl.get(Calendar.MINUTE);
        int second = cl.get(Calendar.SECOND);
        String datetime = String.format("[%04d-%02d-%02dT%02d:%02d:%02d]",year, month, day, hour, minute, second);
        return datetime;
    }
    
    /**
     * Checks the list of previous simulations and returns where in the list a previous simulation with the same parameters occurred.
     * Or a -1 if no such simulation is in the list.
     * @param l The L parameter to find.
     * @param logk The logk to find.
     * @return -1 if there have been no previous simulations with both L and logk, 
     * or the index in {@link #list_data} that corresponds to the given L and logk.
     */
    static int checkList(double l, double e, double logk, int gridsum, double xa, double ya){
        for(int i = 0; i < list_l.size(); i++){
            boolean leq = precisionEquals(l, list_l.get(i), L_DECIMALS);
            boolean e1eq = precisionEquals(e, list_e.get(i), E_DECIMALS);
            boolean e2eq = precisionEquals(logk, list_logk.get(i), LOGK_DECIMALS);
            boolean gseq = list_gridsum.get(i).intValue() == gridsum;
            boolean xaeq = precisionEquals(xa, list_xa.get(i), ADJ_DECIMALS);
            boolean yaeq = precisionEquals(ya, list_ya.get(i), ADJ_DECIMALS);
            if(leq && e1eq && e2eq && gseq && xaeq && yaeq){
                return i;
            }
        }
        return -1;
    }
    
    /**
     * Adds the l, logk, and data to the list of previous simulations.
     * @param l The l parameter.
     * @param logk The logk parameter.
     * @param data The simulated currents.
     */
    static void addToList(double l, double e, double logk, double r, int gridsum, double xa, double ya, double[] data){
        Double dl = l;
        Double de = e;
        Double dlogk = logk;
        Double dr = r;
        Integer gsum = gridsum;
        Double dxa = xa;
        Double dya = ya;
        Double[] ddata = convertToClassDouble(data);
        list_l.add(dl);
        list_e.add(de);
        list_logk.add(dlogk);
        list_gridsum.add(gsum);
        list_data.add(ddata);
        list_r.add(dr);
        list_xa.add(dxa);
        list_ya.add(dya);
    }
    
    /**
     * Method for converting double[] to Double[]
     * @param input The double[] to be converted.
     * @return The equivalent Double[].
     */
    static Double[] convertToClassDouble(double[] input){
        int len = input.length;
        Double[] output = new Double[len];
        for(int i = 0; i < len; i++){
            output[i] = input[i];
        }
        return output;
    }
    
    /**
     * Method for converting Double[] to double[]
     * @param input The Double[] to be converted.
     * @return The equivalent double[].
     */
    static double[] convertToRegularDouble(Double[] input){
        int len = input.length;
        double[] output = new double[len];
        for(int i = 0; i < len; i++){
            output[i] = (double)input[i];
        }
        return output;
    }
    
    /**
     * Determines if d1 and cd2 are equal to one-another to a certain number of decimal points.
     * @param d1 The double to be compared to cd2.
     * @param cd2 The Double to be compared to d1.
     * @param decimals The number of decimal points to which the equality test is done.
     * @return true if round(abs(d1-cd2)*(10**decimals)) == 0, false otherwise.
     */
    static boolean precisionEquals(double d1, Double cd2, int decimals){
        double d2 = (double)cd2;
        double factor = Math.pow(10, decimals);
        double difference = Math.abs(d1 - d2);
        long rounded_difference = Math.round(difference*factor);
        return rounded_difference == 0;
    }
    
    /**
     * Method for finding the best first guess for logk.
     * @param firstl the initial guess for the L parameter.
     * @param verbose flag for determining if extra data logging is desired.
     * @return The first guess for the logarithm of the k parameter.
     * @throws FileNotFoundException 
     */
    static double findFirstLogK(double firstl, double firste, boolean verbose) throws FileNotFoundException, IOException{
        //call the run method to produce a current when the electrode is at distance z and GridData.getCentre() relative to the reactive feature.
		writeReactivityFile(grid, 1.0, 1.0);
        double[] centre = getCentre(1);
        
        double[] kreds = new double[TEST_LOG_K.length];
        double[] koxs = new double[TEST_LOG_K.length];
        
        for(int i = 0; i < TEST_LOG_K.length; i++){
			double test_k0 = Math.pow(10.0, TEST_LOG_K[i]);
            kreds[i] = Math.exp(-ALPHA*F*firste)*test_k0;
            koxs[i] = Math.exp((1-ALPHA)*F*firste)*test_k0;
			
			//impose limits on the rates of reactions to keep them under control.
			
			if(kreds[i] > MAX_KRED){
				double factor = MAX_KRED / kreds[i];
				kreds[i] *= factor;
				koxs[i]  *= factor;
			}
			
			if(koxs[i] > MAX_KOX){
				double factor = MAX_KOX / koxs[i];
				kreds[i] *= factor;
				koxs[i]  *= factor;
			}
        }
        
        Model model = runk(reactivity_mapfile, kreds, koxs, firstl, centre[0], centre[1], 0);
        double[] curr_ox = readData();
        eraseDataFile();
		model = runk(reactivity_mapfile, kreds, koxs, firstl, centre[0], centre[1], 1);
        double[] curr_red = readData();
        eraseDataFile();
        double max_derivative = Math.abs(curr_ox[2] - curr_ox[0]);
        int max_derivative_index = 1;
        for(int i = 2; i < TEST_LOG_K.length -2; i++){
            double derivative_ox = Math.abs(curr_ox[i+1]-curr_ox[i-1]);
			double derivative_red = Math.abs(curr_red[i+1]-curr_red[i-1]);
			double derivative = 0.5*(derivative_ox + derivative_red);
            if(derivative > max_derivative){
                max_derivative = derivative;
                max_derivative_index = i;
            }
        }
        if(verbose){
            writeKFit(KLOGFILE, curr_ox, curr_red);
        }
        return TEST_LOG_K[max_derivative_index];
    }
    
    static double[] powerOfTen(double[] logs){
        double[] output = new double[logs.length];
        for(int i = 0; i < output.length; i++){
            output[i] = Math.pow(10.0, logs[i]);
        }
        return output;
    }
	
    /**
     * Method that controls the lambda-escalation policy for the Levenberg-Marquardt algorithm.
     * @param current_lambda The lambda that was just used.
     * @return The new lambda.
     */
    static double nextLambda(double current_lambda){
        if(current_lambda == 0.0){
            return 1E-4;
        }
        else if(current_lambda == 1E-4){
            return 1E-2;
        }
        else if(current_lambda == 1E-2){
            return 1E-1;
        }
        else if(current_lambda == 1E-1){
            return 0.3;
        }
        else if(current_lambda == 0.3){
            return 1;
        }
        else {
            return 10.0*current_lambda;
        }
    }
    
    /**
     * Computes the sum of squares for the residuals.
     * @param residuals The residuals to be squared and summed.
     * @return The sum of squares.
     */
    static double sumSquare(double[] residuals){
        double sum = 0;
        for(double r : residuals){
            sum += r*r;
        }
        return sum;
    }
    
    /**
     * Rounds value to decimals decimal places.
     * @param value The value to be rounded.
     * @param decimals The number of decimals to round to.
     * @return The rounded value.
     */
    static double round(double value, int decimals){
        return Math.rint(value*Math.pow(10, decimals))/Math.pow(10, decimals);
    }
    
    /**
     * Searches the simulation data for the parameters that have given the lowest sum of square residuals so far.
     * @return The index in the simulation data list that corresponds to the lowest SSR.
     */
    static int findLowestSSR(){
        if(list_l.isEmpty()){
            return -1;
        }
        double[] experimental = true_image;
        int len = list_l.size();
        int low_index = len - 1;
        double[] curr = convertToRegularDouble(list_data.get(low_index));
        double[] residuals = subtract(experimental, curr);
        double low_ssr = sumSquare(residuals);
        for(int i = 0; i < len; i++){
            curr = convertToRegularDouble(list_data.get(i));
            residuals = subtract(experimental, curr);
            double ssr = sumSquare(residuals);
            if(ssr < low_ssr){
                low_index = i;
                low_ssr = ssr;
            }
        }
        return low_index;
    }
    
	/**
     * Apply limit to the l parameter
     * @param l The proposed l parameter.
     * @return The bounded value
     */
    static double applyLLimit(double l){
        double min = 0.1;
        if(l < min){
            l = min;
        }
        return l;
    }
	
    /**
     * Apply limit to the r parameter
     * @param r The proposed r parameter.
     * @return The bounded value
     */
    static double applyRLimit(double r){
        double min = -5;
        double max = 5;
        if(r > max){
            r = max;
        }
        else if(r < min){
            r = min;
        }
        return r;
    }
	
    /**
     * Apply limit to the adjustment parameters
     * @param adj The proposed adjustment parameter.
     * @return The bounded value
     */
    static double applyADJLimit(double adj){
        double min = -5;
        double max = 5;
        if(adj > max){
            adj = max;
        }
        else if(adj < min){
            adj = min;
        }
        return adj;
    }
	
	/**
     * Apply limit to the e parameters
     * @param e The proposed e parameter.
     * @return The bounded value
     */
    static double applyELimit(double e){
        double min = -0.2;
        double max = 0.2;
        if(e > max){
            e = max;
        }
        else if(e < min){
            e = min;
        }
        return e;
    }
	
    /**
     * Apply limit to the logk parameters
     * @param logk The proposed logk parameter.
     * @return The bounded value
     */
    static double applyLOGKLimit(double logk){
        double min = -7;
        double max = -2;
        if(logk > max){
            logk = max;
        }
        else if(logk < min){
            logk = min;
        }
        return logk;
    }
	
    ////////////////////////////////////////////////////////////////////////////
    // CONSTANTS FOR FITTING ///////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    
    /**
     * Exchange rate constant
     */
    static final double K0 = 1;
    
    /*
    Lists for storing simulation history (to avoid redundant simulations)
    */
    
    /**
     * List for holding the l-parameters that have been simulated
     */
    static LinkedList<Double> list_l;
    
    /**
     * List for holding the logk parameters that have been simulated
     */
    static LinkedList<Double> list_e;
    
    /**
     * List for holding the logk parameters that have been simulated
     */
    static LinkedList<Double> list_logk;
    
    /**
     * List for holding the total number of grid points that are switched 'on'.
     * Since the possible grid configurations are controlled by erosion or dilation of the initial grid configuration,
     * The pixels will turn 'on' or 'off' in a very specific order.
     */
    static LinkedList<Integer> list_gridsum;
    
    /**
     * List for holding r-parameters that have been simulated
     */
    static LinkedList<Double> list_r;
    
    /**
     * List for holding the current results for each secm image simulation.
     */
    static LinkedList<Double[]> list_data;
    
    /**
    * List for holding x-adjustment-parameters that have been simulated
    */
    static LinkedList<Double> list_xa;
	
    /**
    * List for holding y-adjustment-parameters that have been simulated
    */
    static LinkedList<Double> list_ya;
	
    /*
    Fitting status symbols
    */
    
    /**
     * Indicates that the fitting ended in convergence.
     */
    static final int EXECUTED_OK = 0;
    
    /**
     * Indicates that the fitting ended due to reaching the maximum number of iterations
     */
    static final int MAX_ITERATIONS_REACHED = 1;
	
	static final double MAX_KOX = 1.0;
	
	static final double MAX_KRED = 1.0;
    
    /**
     * Indicates that the fitting method ended due to the maximum value of lambda being reached prior to convergence.
     */
    static final int MAX_LAMBDA_REACHED = 2;
    
    
    /*
    Grid-handling methods and fields
    */
    
    /**
     * X-size of pixel {@link #grid} must be odd
     */
    static final int XSIZE = 201;
    
    /**
     * Y-size of pixel {@link #grid} must be odd
     */
    static final int YSIZE = 201;
    
    /**
     * Obtains the most central sampling point over the reactive feature for obtaining the initial logk guess.
     * This is achieved by scoring each point by the sum of the reciprocal distance squared to each of the 'on' pixels.
     * @return The position of the most central sampling point as an array {x, y}.
     */
    static double[] getCentre(int selector){
        
        double bestx = 0;
        double besty = 0;
        double bestScore = 0;
        
        for(int sx = 0; sx < sample_xs.length; sx++){
            int ix = sample_xs[sx] - min_x;
            for(int sy = 0; sy < sample_ys.length; sy++){
                int iy = sample_ys[sy] - min_y;
                double score = 0;
                for(int x = 0; x < grid.length; x++){
                    for(int y = 0; y < grid[0].length; y++){
                        if((x != ix || y != iy) && grid[x][y] == selector){
                            double distsq = (ix-x)*(ix-x) + (iy-y)*(iy-y);
                            score += 1.0/distsq;
                        }
                    }
                }
                if(score > bestScore){
                    bestx = physical_xs[sx];
                    besty = physical_ys[sy];
                    bestScore = score;
                }
            }
        }
        return new double[]{bestx, besty};
    }
    
    /**
     * The experimental SECM currents
     */
    static double[] true_image;
    
    /**
     * The experimental x-coordinates in meters
     */
    static double[] physical_xs;
    
    /**
     * The experimental y-coordinates in meters
     */
    static double[] physical_ys;
    
	/**
	 * Switch that determines the ME behaviour. (0 for oxidising, 1 for reducing)
	 */
	static double[] me_o0r1;
	
    /**
     * The grid_switches x-indexes
     */
    static int[] sample_xs;
    
    /**
     * The grid_switches y-indexes
     */
    static int[] sample_ys;
    
    /**
     * The pixel grid_switches
     */
    static int[][] grid;
    
    /**
     * the minimum grid_switches x-index
     */
    static int min_x;
    
    /**
     * the minimum grid_switches y-index
     */
    static int min_y;
    
    /**
     * The filepath for the instruction file.
     */
    static String true_file;
    
    /*
    IO-handling methods and fields
    */
    
    /**
     * Reads-in the current data from the data file that is produced by the simulation.
     * @return The currents at each point in the same order as {@link #sample_xs}, {@link #sample_ys} and {@link #true_image}.
     * COMSOL prefixes comment lines with '%', so lines starting with '%' in the file are ignored.
     * @throws FileNotFoundException 
     */
    static double[] readData() throws FileNotFoundException{
        File f = new File("data.txt");
        Scanner s = new Scanner(f);
        int count = 0;
        while(s.hasNextLine()){
            String ln = s.nextLine();
            if(!ln.startsWith("%")){
                count ++;
            }
        }
        s.close();
        double[] data = new double[count];
        s = new Scanner(f);
        count = 0;
        while(s.hasNextLine()){
            String ln = s.nextLine();
            if(!ln.startsWith("%")){
                String[] linedata = ln.trim().split("\\s+");//split up the columns
                data[count] = Double.parseDouble(linedata[linedata.length - 1]);//take the data from the righmost column
                count ++;
            }
        }
        return data;
    }
    
    /**
     * Erases the data file produced by the simulation, writing a '%' to the file.
     * @throws FileNotFoundException 
     */
    static void eraseDataFile() throws FileNotFoundException{
        File f = new File("data.txt");
        PrintWriter pw = new PrintWriter(f);
            pw.print("%");
        pw.close();
    }
    
    /**
     * Logs iteration details into {@link #LOGFILE}.
     * @param iteration_num The number of the iteration.
     * @param diagonal The diagonal of the DTD matrix
     * @param lambda The lambda value
     * @param labels The labels for the parameters
     * @param params The values of the parameters
     * @param ssr The sum of square residuals for the parameters.
     * @param accepted Whether or not the parameters had an acceptable SSR.
     * @throws FileNotFoundException
     * @throws IOException 
     */
    static void logIteration(int iteration_num, double[] diagonal, double lambda, String[] labels, double[] params, double ssr, boolean accepted) throws FileNotFoundException, IOException{
        File f = new File(LOGFILE);
        PrintWriter pw = new PrintWriter(new FileWriter(f, true));
            pw.append("\nIteration: " + iteration_num);
            pw.append("\n\t Diagonal: " + diagonal[0]);
            for(int i = 1; i < labels.length; i++){
                pw.append("," + diagonal[i]);
            }
        pw.close();
        logLambda(lambda, labels, params, ssr, accepted);
    }
    
    /**
     * Logs details regarding an attempt with a new lambda to {@link #LOGFILE}.
     * @param lambda The lambda value being used.
     * @param labels The parameter labels.
     * @param params The parameter values.
     * @param ssr The sum of square residuals.
     * @param accepted Whether or not the lambda attempt had an acceptable SSR.
     * @throws FileNotFoundException
     * @throws IOException 
     */
    static void logLambda(double lambda, String[] labels, double[] params, double ssr, boolean accepted) throws FileNotFoundException, IOException{
        File f = new File(LOGFILE);
        PrintWriter pw = new PrintWriter(new FileWriter(f, true));
            if(lambda >= 0){
				pw.append("\n\tLAMBDA: " + lambda);
			}
			else{
				pw.append("\n\tFound previous simulation with lower SSR");
			}
            for(int i = 0; i < labels.length; i++){
                pw.append("\n\t\t" + labels[i] + ": " + params[i]);
            }
            pw.append("\n\tSum. Square Residuals: " + ssr);
            if(accepted){
                pw.append(" (Accepted)");
            }
            else{
                pw.append(" (Rejected)");
            }
        pw.close();
    }
    
    /**
     * Logs the SSR for the initial parameters to {@link #LOGFILE}.
     * @param fname The filename containing the experimental data that is being fit. 
     * (Note that this parameter is not for the name of the file that this logging information is being written into)
     * @param labels The labels for the parameters.
     * @param params The values of the parameters.
     * @param ssr The sum of square residuals when using the parameters.
     * @throws FileNotFoundException
     * @throws IOException 
     */
    static void logInitialGuesses(String fname, String[] labels, double[] params, double ssr) throws FileNotFoundException, IOException{
        File f = new File(LOGFILE);
        PrintWriter pw = new PrintWriter(f);
            pw.append("\nInitial guesses:");
            pw.append("\n\tFile: " + fname);
            for(int i = 0; i < labels.length; i++){
                pw.append("\n\t\t" + labels[i] + ": " + params[i]);
            }
            pw.append("\n\tSum. Square Residuals: " + ssr);
        pw.close();
    }
    
    /**
     * Logs the end result for the fitting to {@link #LOGFILE}.
     * @param condition The reason for the fit ending.
     * @throws FileNotFoundException
     * @throws IOException 
     * @see #EXECUTED_OK
     * @see #MAX_ITERATIONS_REACHED
     * @see #MAX_LAMBDA_REACHED
     */
    static void logEndCondition(int condition) throws FileNotFoundException, IOException{
        File f = new File(LOGFILE);
        PrintWriter pw = new PrintWriter(new FileWriter(f, true));
            switch(condition){
                case EXECUTED_OK:
                    pw.append("\nPROCESS CONVERGED.");
					break;
                case MAX_ITERATIONS_REACHED:
                    pw.append("\nPROCESS STOPPED PREMATURELY AFTER " + MAX_ITERATIONS + " ITERATIONS.");
					break;
                case MAX_LAMBDA_REACHED:
                    pw.append("\nPROCESS STOPPED DUE TO MAXIMUM LAMBDA BEING REACHED.");
            }
        pw.close();
    }
    
    /**
     * Produces a copy of the file at original_filepath, replacing the currents at the sampled points with the currents in {@link #true_image}.
     * This copy will be found at new_filepath.
     * This method is intended to be used to easily generate simulated secm images for testing purposes.
	 * THIS METHOD ASSUMES true_image is separated by me_o0r1 before being separated by x and y. i.e. true image is the oxidized or reduced image appended after the reduced or oxidised image
     * @param original_filepath The path to a correctly formatted experimental control file.
     * @param new_filepath The path to the file that is to be created/written-to.
     * @throws FileNotFoundException
     * @throws IOException 
     */
    static void writeSECMInfo(String original_filepath, String new_filepath) throws FileNotFoundException, IOException{
        File originalfile = new File(original_filepath);
        File newfile = new File(new_filepath);
		int half = true_image.length / 2;
        if(original_filepath.equals(new_filepath)){
            throw new IOException("Cannot read and write to the same file");
        }
        Scanner s = new Scanner(new BufferedReader(new FileReader(originalfile)));
        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(newfile)));
        int index = 0;
        String sep = ",";
        
        if(s.hasNextLine()){
            String fl = s.nextLine();
            pw.print(fl);
            if(fl.equalsIgnoreCase("##ENCODING: csv")){
                sep = ",";
            }
            else if(fl.equalsIgnoreCase("##ENCODING: tsv")){
                sep = "\t";
            }
            else{
                throw new FileNotFoundException("File incorrectly formatted.");
            }
        }
        else{
            throw new FileNotFoundException("File incorrectly formatted.");
        }
        
        while(s.hasNextLine()){
            String line = s.nextLine();
            if(!line.startsWith("#")){
                String[] linesplit = line.split(sep);
                int x = Integer.parseInt(linesplit[0]);
                int y = Integer.parseInt(linesplit[1]);
                if(x == sample_xs[index] && y == sample_ys[index]){
                    //copy-over all of the line save the current
					if(me_o0r1[index] < 0.1){
						pw.print("\n" + linesplit[0] + sep + linesplit[1] + sep + linesplit[2] + sep + linesplit[3] + sep + linesplit[4] + sep + true_image[index] + sep + true_image[index + half]);
					}
					else{
						pw.print("\n" + linesplit[0] + sep + linesplit[1] + sep + linesplit[2] + sep + linesplit[3] + sep + linesplit[4] + sep + true_image[index + half] + sep + true_image[index]);
					}
                    index ++;
                    if(index >= sample_xs.length){
                        index = sample_xs.length - 1;
                    }
                }
                else{
                    pw.print("\n" + line);//copy-over non-sampled lines
                }
            }
            else{
                pw.print("\n" + line);//copy-over commented lines
            }
        }
        
        s.close();
        pw.close();
    }
    
    /**
     * Reads through the control file, pulling out information for the reactivity grid, the sampled points and the experimentally observed secm currents.
     * This method instantiates:
     * <ul>
     * <li>{@link #grid}</li>
     * <li>{@link #sample_xs}</li>
     * <li>{@link #sample_ys}</li>
     * <li>{@link #physical_xs}</li>
     * <li>{@link #physical_ys}</li>
     * <li>{@link #min_x}</li>
     * <li>{@link #min_y}</li>
     * </ul>
     * @param filepath The path to the control file.
     * @throws FileNotFoundException 
     */
    static void readSECMInfo(String filepath) throws FileNotFoundException{
        File f = new File(filepath);
        true_file = filepath;
        int minx = 0;
        int miny = 0;
        int max_x = 0;
        int max_y = 0;
        double[] trueimage;
        double[] physicalxs;
        double[] physicalys;
		double[] localme_o0r1;
        int[] samplexs;
        int[] sampleys;
        String sep = ",";
        
        Scanner s = new Scanner(new BufferedReader(new FileReader(f)));
            if(s.hasNextLine()){
                String fl = s.nextLine();
                if(fl.equalsIgnoreCase("##ENCODING: csv")){
                    sep = ",";
                }
                else if(fl.equalsIgnoreCase("##ENCODING: tsv")){
                    sep = "\t";
                }
                else{
                    throw new FileNotFoundException("File incorrectly formatted.");
                }
            }
            else{
                throw new FileNotFoundException("File incorrectly formatted.");
            }
            
            //X header
            int xstart = 0;
            int xstep = 0;
            int xnum = 0;
            if(s.hasNextLine()){
                String fl = s.nextLine();
            }
            else{
                throw new FileNotFoundException("File incorrectly formatted.");
            }
            if(s.hasNextLine()){
                String fl = s.nextLine();
                String[] tokens = fl.substring(1).trim().split(",");
                xstart = Integer.parseInt(tokens[0]);
                xstep = Integer.parseInt(tokens[1]);
                xnum = Integer.parseInt(tokens[2]);
            }
            else{
                throw new FileNotFoundException("File incorrectly formatted.");
            }
            
            //Y header
            int ystart = 0;
            int ystep = 0;
            int ynum = 0;
            if(s.hasNextLine()){
                String fl = s.nextLine();
            }
            else{
                throw new FileNotFoundException("File incorrectly formatted.");
            }
            if(s.hasNextLine()){
                String fl = s.nextLine();
                String[] tokens = fl.substring(1).trim().split(",");
                ystart = Integer.parseInt(tokens[0]);
                ystep = Integer.parseInt(tokens[1]);
                ynum = Integer.parseInt(tokens[2]);
            }
            else{
                throw new FileNotFoundException("File incorrectly formatted.");
            }
            
            //read through the main data file to pull out trueimage, physicalxs, physicalys, samplexs, sampleys, minx, miny
            minx = xstart;
            max_x = xstart;
            miny = ystart;
            max_y = ystart;
            
            int asize = (xnum)*(ynum);
            int present_index = 0;
            
            trueimage = new double[asize*2];
            physicalxs = new double[asize*2];
            physicalys = new double[asize*2];
			localme_o0r1 = new double[asize*2];
            samplexs = new int[asize*2];
            sampleys = new int[asize*2];
            
            while(s.hasNextLine()){
                String line = s.nextLine();
                if(!line.startsWith("#")){
                    String[] linesplit = line.split(sep);
                    int x = Integer.parseInt(linesplit[0]);
                    int y = Integer.parseInt(linesplit[1]);
                    double px = Double.parseDouble(linesplit[3]);
                    double py = Double.parseDouble(linesplit[4]);
                    double curr_mox = Double.parseDouble(linesplit[5]);//current when ME is oxidising
					double curr_mre = Double.parseDouble(linesplit[6]);//current when ME is reducing
                    
                    if(x < minx){
                        minx = x;
                    }
                    else if(x > max_x){
                        max_x = x;
                    }
                    
                    if(y < miny){
                        miny = y;
                    }
                    else if(y > max_y){
                        max_y = y;
                    }
                    
                    //check if x,y is one of the points to be sampled
                    boolean xvalid = (x >= xstart) && (x < xstart + xnum*xstep) && ((x-xstart)%xstep == 0);
                    boolean yvalid = (y >= ystart) && (y < ystart + ynum*ystep) && ((y-ystart)%ystep == 0);
                    if(xvalid && yvalid && present_index < asize){
                        trueimage[present_index] = curr_mox;
						trueimage[present_index + asize] = curr_mre;
						localme_o0r1[present_index] = 0;
						localme_o0r1[present_index + asize] = 1;
                        physicalxs[present_index] = px;
						physicalxs[present_index + asize] = px;
                        physicalys[present_index] = py;
						physicalys[present_index + asize] = py;
                        samplexs[present_index] = x;
						samplexs[present_index + asize] = x;
                        sampleys[present_index] = y;
						sampleys[present_index + asize] = y;
                        present_index ++;
                    }
                    
                }
            }
        s.close();
	    if(present_index != asize){
		throw new FileNotFoundException("Expected " + asize + " points. Found: " + present_index + ".");
	    }
        //second read of the file to get grid_switches[][]
        int[][] grid_switches = new int[max_x - minx + 1][max_y - miny + 1];
        s = new Scanner(new BufferedReader(new FileReader(f)));
            while(s.hasNextLine()){
                String line = s.nextLine();
                if(!line.startsWith("#")){
                    String[] linesplit = line.split(sep);
                    int x = Integer.parseInt(linesplit[0]);
                    int y = Integer.parseInt(linesplit[1]);
                    int rs = Integer.parseInt(linesplit[2]);
                    
                    grid_switches[x-minx][y-miny] = rs;
                    
                }
            }
        s.close();
        
        true_image = trueimage;
        physical_xs = physicalxs;
        physical_ys = physicalys;
        sample_xs = samplexs;
        sample_ys = sampleys;
        min_x = minx;
        min_y = miny;
        grid = grid_switches;
		me_o0r1 = localme_o0r1;
        
    }
    
    /**
     * Writes the current, derivatives and residuals to a file for plotting purposes.
     * @param fname The name of the file to be written to.
     * @param x The x-positions of the sampled points.
     * @param y The y-positions of the sampled points.
     * @param current The currents.
     * @param dl The Partial derivatives of current with respect to L.
     * @param dlogk The partial derivatives of current with respect to logk.
     * @param residual The residuals.
     * @throws FileNotFoundException 
     */
    static void writeIteration(String fname, double[] x, double[] y, double[] me_o0_r1, double[] current, double[] dl, double[] de, double[] dlogk, double[] dr, double[] dxa, double[] dya, double[] residual) throws FileNotFoundException{
        File f = new File(fname);
        PrintWriter pw = new PrintWriter(f);
            pw.print("#x [m], y [m], ox:0 | red:1, i [A], dL [A], de [A/V], dlogk [A], dr [A], dxa [A], dya [A], residual [A]");
            for(int i = 0; i < x.length; i ++){
                pw.print("\n" + x[i] + "," + y[i] + "," + me_o0_r1[i] + "," + current[i] + "," + dl[i] + "," + de[i] + "," + dlogk[i] + "," + dr[i] + "," + dxa[i] + "," + dya[i] + "," + residual[i]);
            }
        pw.close();
    }

    /**
     * Writes the logk-current curve data to a file.
     * @param fname The file to be written to.
     * @param current The currents from the logk simulations.
     * @throws FileNotFoundException 
     */
    static void writeKFit(String fname, double[] current_ox, double[] current_red) throws FileNotFoundException{
        File f = new File(fname);
        PrintWriter pw = new PrintWriter(f);
            pw.print("#log10(k/1[m/s]), i_ox [A], i_red [A]");
            for(int i = 0; i < current_ox.length; i ++){
                pw.print("\n" + TEST_LOG_K[i] + "," + current_ox[i] + "," + current_red[i]);
            }
        pw.close();
    }
    
    /**
     * Writes the reactivity information to {@link #reactivity_mapfile}
     * @param de_grid The map of pixels that are switched on and off
     * @param reactivity The rate constant of the activated pixels in [m/s]
     * @throws IOException
     * @throws FileNotFoundException 
     */
    static void writeReactivityFile(int[][] de_grid, double kred, double kox) throws IOException, FileNotFoundException{
        double[][] kredmap = new double[grid.length][grid[0].length];
        double[][] koxmap = new double[grid.length][grid[0].length];
        for(int x = 0; x < grid.length; x++){
            for(int y = 0; y < grid[0].length; y++){
                int pv = de_grid[x][y];
                switch (pv) {
                    case 1:
                        kredmap[x][y] = kred;
                        koxmap[x][y] = kox;
                        break;
//                    case 2:
//                        kredmap[x][y] = kred2;
//                        koxmap[x][y] = kox2;
//                        break;
                    default:
                        kredmap[x][y] = 0.0;
                        koxmap[x][y] = 0.0;
                        break;
                }
            }
        }
        File originalfile = new File(true_file);
        File newfile = new File(reactivity_mapfile);
        if(true_file.equals(reactivity_mapfile)){
            throw new IOException("Cannot read and write to the same file");
        }
        newfile.createNewFile();
        Scanner s = new Scanner(new BufferedReader(new FileReader(originalfile)));
        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(newfile)));
        String sep = ",";
        
        if(s.hasNextLine()){
            String fl = s.nextLine();
            if(fl.equalsIgnoreCase("##ENCODING: csv")){
                sep = ",";
            }
            else if(fl.equalsIgnoreCase("##ENCODING: tsv")){
                sep = "\t";
            }
            else{
                throw new FileNotFoundException("File incorrectly formatted.");
            }
        }
        else{
            throw new FileNotFoundException("File incorrectly formatted.");
        }
        
        boolean first = true;
        while(s.hasNextLine()){
            String line = s.nextLine();
            if(!line.startsWith("#")){
                String[] linesplit = line.split(sep);
                int x = Integer.parseInt(linesplit[0]);
                int y = Integer.parseInt(linesplit[1]);
                //copy-over the relevant parts of the file
                if(!first){
                    pw.print("\n" + linesplit[3] + "," + linesplit[4] + "," + kredmap[x][y] + "," + koxmap[x][y]);
                }
                else{
                    pw.print(linesplit[3] + "," + linesplit[4] + "," + kredmap[x][y] + "," + koxmap[x][y]);
                    first = false;
                }
            }
        }
        
        s.close();
        pw.close();
    }
    
    /**
     * Fetches the current working directory.
     * @return The absolute file path of the directory from which this program is executing, ".".
     */
    static String getCWD(){
		File rmf = new File(reactivity_mapfile);
		File cwd = rmf.getAbsoluteFile().getParentFile();
		return cwd.getPath();
    }
    
    /**
     * Converts an array of doubles to a space separated String.
     * @param a The array of doubles to be converted.
     * @return A space separated String containing all of the elements of a.
     */
    static String toString(double[] a){
        if(a.length > 1){
            String sep = " ";
            String out = "" + a[0];
            for(int i = 1; i < a.length; i++){
                out = out + sep + a[i];
            }
            return out;
        }
        else if(a.length == 1){
            return "" + a[0];
        }
        else{
            return "";
        }
    }
    
    /**
     * The file to which logging information will be written.
     */
    static final String LOGFILE = "fit.log";
    
    /**
     * The file to which logk-current curves will be written.
     */
    static final String KLOGFILE = "k-curve.csv";
    
    /**
     * The file to which the reactivity map will be written and read from.
     */
    static String reactivity_mapfile = "func.csv";
    
    /*
    Linear algebra-handling methods and fields
     */
    /**
     * Inverts the given matrix
     * @param a the matrix to be inverted
     * @return the inverse of a, a<sup>-1</sup>. a<sup>-1</sup>*a = a*a<sup>-1</sup>=I
     */
    static double[][] invert(double[][] a) throws NumberFormatException{
        double det = determinant(a);
        if(!Double.isFinite(det) || det == 0.0){
            throw new NumberFormatException("Matrix a is singular.");
        }
        else{
            double[][] cofactors = new double[a.length][a[0].length];
            for(int r = 0; r < a.length; r ++){
                for(int c = 0; c < a.length; c ++){
                    cofactors[r][c] = cofactor(a, r, c);
                }
            }
            return multiply(transpose(cofactors), 1.0/det);
        }
    }
    
    /**
     * Multiplies matrices a*b
     * @param a The left matrix. a[r][c] (r = row number; c = column number)
     * @param b The right matrix. b[r][c] (r = row number; c = column number)
     * @return The product matrix. p[r][c] (r = row number; c = column number)
     * @throws NumberFormatException 
     */
    static double[][] multiply(double[][] a, double[][] b) throws NumberFormatException{
        double[][] prod = new double[a.length][b[0].length];
        if(a[0].length != b.length){
            throw new NumberFormatException("columnspace of a, " + a[0].length + ", and rowspace of b, " + b.length + ", must be the same size.");
        }
        for(int r = 0; r < a.length; r++){
            for(int c = 0; c < b[0].length; c++){
                prod[r][c] = 0;
                for(int i = 0; i < b.length; i++){
                    prod[r][c] += a[r][i]*b[i][c];
                }
            }
        }
        return prod;
    }
    
    /**
     * Operates the matrix a on the vector v. a*v
     * @param a The left matrix. a[r][c] (r = row number; c = column number)
     * @param v The right vector. v[r] (r = row number)
     * @return The resulting product vector. p[r] (r = row number)
     * @throws NumberFormatException 
     */
    static double[] multiply(double[][]a, double[] v) throws NumberFormatException{
        if(a[0].length != v.length){
            throw new NumberFormatException("columnspace of a, " + a[0].length + ", and rowspace of v, " + v.length + ", must be the same size.");
        }
        double[] prod = new double[a.length];
        for(int r = 0; r < a.length; r++){
            prod[r] = 0;
            for(int c = 0; c < v.length; c++){
                prod[r] += a[r][c]*v[c];
            }
        }
        return prod;
    }
    
    /**
     * Multiplies a matrix a with a scalar s.
     * @param a The a matrix
     * @param s The s scalar
     * @return The product matrix
     */
    static double[][] multiply(double[][] a, double s){
        double[][] mul = new double[a.length][a[0].length];
        for(int r = 0; r < a.length; r++){
            for(int c = 0; c < a[0].length; c++){
                mul[r][c] = s*a[r][c];
            }
        }
        return mul;
    }
    
    /**
     * Multiplies a vector v with a scalar s.
     * @param v The a vector
     * @param s The s scalar
     * @return The product vector
     */
    static double[] multiply(double[] v, double s){
        double[] mul = new double[v.length];
        for(int i = 0; i < mul.length; i++){
            mul[i] = v[i]*s;
        }
        return mul;
    }
    
    /**
     * Adds two matrices a+b
     * @param a The left matrix. a[r][c] (r = row number; c = column number)
     * @param b The right matrix. b[r][c] (r = row number; c = column number)
     * @return The sum matrix. s[r][c] (r = row number; c = column number)
     * @throws NumberFormatException 
     */
    static double[][] add(double[][] a, double[][] b) throws NumberFormatException{
        if(a.length != b.length || a[0].length != b[0].length){
            throw new NumberFormatException("a and b must be the same size!");
        }
        double[][] sum = new double[a.length][a[0].length];
        for(int r = 0; r < a.length; r++){
            for(int c = 0; c < a[0].length; c++){
                sum[r][c] = a[r][c] + b[r][c];
            }
        }
        return sum;
    }
    
    /**
     * Subtracts two matrices a-b
     * @param a The left matrix. a[r][c] (r = row number; c = column number)
     * @param b The right matrix. b[r][c] (r = row number; c = column number)
     * @return The difference matrix. d[r][c] (r = row number; c = column number)
     * @throws NumberFormatException 
     */
    static double[][] subtract(double[][] a, double[][] b) throws NumberFormatException{
        if(a.length != b.length || a[0].length != b[0].length){
            throw new NumberFormatException("a and b must be the same size!");
        }
        double[][] dif = new double[a.length][a[0].length];
        for(int r = 0; r < a.length; r++){
            for(int c = 0; c < a[0].length; c++){
                dif[r][c] = a[r][c] - b[r][c];
            }
        }
        return dif;
    }
    
    /**
     * Adds two vectors a+b
     * @param a The left vector. a[i]
     * @param b The right vector. b[i]
     * @return The sum vector. s[i]
     * @throws NumberFormatException 
     */
    static double[] add(double[] a, double[] b) throws NumberFormatException{
        if(a.length != b.length){
            throw new NumberFormatException("a and b must be the same size!");
        }
        double[] sum = new double[a.length];
        for(int r = 0; r < a.length; r++){
            sum[r] = a[r] + b[r];
        }
        return sum;
    }
    
    /**
     * Subtracts two vectors a-b
     * @param a The left vector. a[i]
     * @param b The right vector. b[i]
     * @return The difference vector. d[i]
     * @throws NumberFormatException 
     */
    static double[] subtract(double[] a, double[] b) throws NumberFormatException{
        if(a.length != b.length){
            throw new NumberFormatException("a and b must be the same size!");
        }
        double[] dif = new double[a.length];
        for(int r = 0; r < a.length; r++){
            dif[r] = a[r] - b[r];
        }
        return dif;
    }
    
    /**
     * Swaps the column and rowspace of a.
     * @param a The matrix. a[r][c] (r = row number; c = column number) 
     * @return The transposed matrix. t[r][c] (r = row number; c = column number) 
     */
    static double[][] transpose(double[][] a){
        double[][] transp = new double[a[0].length][a.length];
        for(int r = 0; r < a.length; r ++){
            for(int c = 0; c < a[0].length; c ++){
                transp[c][r] = a[r][c];
            }
        }
        return transp;
    }
    
    /**
     * appends v to the right of a
     * @param a the matrix to which v is to be appended
     * @param v the vector to append
     * @return a copy of a with column v added to the right
     * @throws NumberFormatException 
     */
    static double[][] appendColumn(double[][] a, double[] v) throws NumberFormatException{
        if(a.length != v.length){
            throw new NumberFormatException("Rowspace of a, " + a.length + ", and v, " + v.length + ", must be the same");
        }
        double[][] augment = new double[a.length][a[0].length + 1];
        for(int r = 0; r < a.length; r++){
            for(int c = 0; c < a[0].length; c++){
                augment[r][c] = a[r][c];
            }
            augment[r][a[0].length] = v[r];
        }
        return augment;
    }
    
    /**
     * appends v to the right of a
     * @param a the vector to which v is to be appended
     * @param v the vector to append
     * @return a copy of a with column v added to the right
     * @throws NumberFormatException 
     */
    static double[][] appendColumn(double[] a, double[] v) throws NumberFormatException{
        if(a.length != v.length){
            throw new NumberFormatException("Rowspace of a, " + a.length + ", and v, " + v.length + ", must be the same");
        }
        double[][] augment = new double[a.length][2];
        for(int r = 0; r < a.length; r++){
            augment[r][0] = a[r];
            augment[r][1] = v[r];
        }
        return augment;
    }
    
    /**
     * Returns the nxn identity matrix, I.
     * @param n the dimension of the identity matrix
     * @return the nxn identity matrix
     */
    static double[][] identity(int n){
        double[][] ident = new double[n][n];
        for(int r = 0; r < n; r++){
            for(int c = 0; c < n; c++){
                if(r==c){
                    ident[r][c] = 1;
                }
                else{
                    ident[r][c] = 0;
                }
            }
        }
        return ident;
    }
    
    /**
     * Computes the determinant of a matrix, a.
     * @param a The matrix a[r][c] (r = row number; c = column number)
     * @return The determinant of a.
     * @throws NumberFormatException 
     * @see #inner_determinant(double[][]) 
     */
    static double determinant(double[][] a) throws NumberFormatException{
        if(a.length == a[0].length){
            return inner_determinant(a);
        }
        else{
            throw new NumberFormatException("Determinants may only be taken of square matrices.");
        }
    }
    
    /**
     * Please use this method's wrapper, {@link #determinant(double[][])} as this function doe not perform any input checking.
     * @param a
     * @return 
     */
    static double inner_determinant(double[][] a){
        if(a.length == 2){//2x2 matrix
            return a[0][0]*a[1][1] - a[0][1]*a[1][0];
        }
        else if(a.length == 1){//1x1 matrix
            return a[0][0];
        }
        else{
            double sum = 0;
            for(int c = 0; c < a[0].length; c++){
                if(c%2 == 0){
                    sum += a[0][c]*inner_determinant(minor(a,0,c));
                }
                else{
                    sum -= a[0][c]*inner_determinant(minor(a,0,c));
                }
            }
            return sum;
        }
    }
    
    /**
     * Generates the cofactor for the element of matrix a at r,c
     * @param a The matrix a[r][c] (r = row number; c = column number)
     * @param r The row index
     * @param c The column index
     * @return The cofactor for constructing the matrix of cofactors for a.
     * @throws NumberFormatException 
     */
    static double cofactor(double[][] a, int r, int c) throws NumberFormatException{
        double sign = -1.0;
        if ((r+c)%2 == 0){
            sign = 1.0;
        }
        return sign*determinant(minor(a, r, c));
    }
    
    /**
     * Returns a copy of a that omits row r and column c.
     * @param a the original matrix
     * @param r the row to be omitted
     * @param c the column to be omitted
     * @return
     */
    static double[][] minor(double[][] a, int r, int c){
        double[][] red = new double[a.length - 1][a[0].length - 1];
        
        for(int ri = 0; ri < a.length - 1; ri++){
            int ru = ri;//the r to be used
            if(ri >= r){
                ru ++;//if we are at or beyond the row to be omitted, increment the r to be read from a
            }
            for(int ci = 0; ci < a[0].length - 1; ci++){
                int cu = ci;//the c to be used
                if(ci >= c){
                    cu ++;//if we are at or beyond the column to be omitted, increment the c to be read from a
                }
                red[ri][ci] = a[ru][cu];
            }
        }
        
        return red;
    }
    
    /* 
    Erosion and dilation
    */
    /**
     * Applies a dilation or erosion by |amount|
     * If amount is negative, an erosion will be carried-out, otherwise, a dilation will.
     * @param amount the amount by which the grid is to be dilated or eroded in grid-space (i.e. the distance between adjacent grid elements is 1.
     */
    static int[][] applyDilationErosion(double amount){
        if(amount < 0){
            return erodeGrid(-amount);
        }
        else{
            return dilateGrid(amount);
        }
    }
    /**
     * dilates grid by amount, storing the result to edited_grid
     * @param amount the amount by which the grid is to be dilated in grid-space (i.e. the distance between adjacent grid elements is 1.
     */
    static int[][] dilateGrid(double amount){
        int upperbound = (int)Math.ceil(amount);
        double amountsq = (amount + 0.5) * (amount + 0.5);// adding 0.5 due to the dilation being from the center of the pixel
        //(re)-initialize the edited_grid
        int[][] edited_grid = new int[grid.length][grid[0].length];
        for(int x = 0; x < grid.length; x++){
            for(int y = 0; y < grid[0].length; y++){
                edited_grid[x][y] = grid[x][y];
            }
        }
        if(amount >= 0.2){ // don't bother with the dilation if the amount is too small to change anything
            //iterate through every point
            for(int x = 0; x < grid.length; x++){
                for(int y = 0; y < grid[0].length; y++){
                    if(grid[x][y] == 0){
                        //if the value of the grid point can be changed, iterate through its neighboring points
                        int[][] subdivisions1 = new int[5][5];
                        int[][] subdivisions2 = new int[5][5];
                        for(int xx = - upperbound; xx <= upperbound; xx++){
                            for(int yy = - upperbound; yy <= upperbound; yy++){
                                if(getGrid(x+xx, y+yy) == 1){
                                    /*if the point at (x+xx,y+yy) is different from 
                                    the one at (x,y), subdivide the point at (x,y)
                                    into a 5x5 grid and see how many of these grid
                                    points are within the erosion or dilation distance
                                    */
                                    for(int xxx = 0; xxx < 5; xxx ++){
                                        for(int yyy = 0; yyy < 5; yyy ++){
                                            double x_subdist = ((double)xxx)*0.2 - 0.4;
                                            double y_subdist = ((double)yyy)*0.2 - 0.4;
                                            double xdist = ((double)xx) + x_subdist;
                                            double ydist = ((double)yy) + y_subdist;
                                            double distancesq = xdist*xdist + ydist*ydist;
                                            if(distancesq < amountsq){
                                                subdivisions1[xxx][yyy] = 1;
                                            }
                                        }
                                    }
                                }
                                else if(getGrid(x+xx, y+yy) == 2){
                                    /*if the point at (x+xx,y+yy) is different from 
                                    the one at (x,y), subdivide the point at (x,y)
                                    into a 5x5 grid and see how many of these grid
                                    points are within the erosion or dilation distance
                                    */
                                    for(int xxx = 0; xxx < 5; xxx ++){
                                        for(int yyy = 0; yyy < 5; yyy ++){
                                            double x_subdist = ((double)xxx)*0.2 - 0.4;
                                            double y_subdist = ((double)yyy)*0.2 - 0.4;
                                            double xdist = ((double)xx) + x_subdist;
                                            double ydist = ((double)yy) + y_subdist;
                                            double distancesq = xdist*xdist + ydist*ydist;
                                            if(distancesq < amountsq){
                                                subdivisions2[xxx][yyy] = 1;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        /*compute the number of subdivisions at (x,y) that should change.
                        If they are in the majority, flip the pixel.
                        */
                        int sum1 = 0;
                        int sum2 = 0;
                        for(int i = 0; i < 5; i++){
                            for(int ii = 0; ii < 5; ii++){
                                sum1 += subdivisions1[i][ii];
                                sum2 += subdivisions2[i][ii];
                            }
                        }
                        if(sum1 > 12){
                            edited_grid[x][y] = 1;
                        }
                        else if(sum2 > 12){
                            edited_grid[x][y] = 2;
                        }
                    }
                }
            }
        }
        return edited_grid;
    }
    
    /**
     * erodes grid by amount, storing the result to edited_grid
     * @param amount the amount by which the grid is to be eroded in grid-space (i.e. the distance between adjacent grid elements is 1.
     */
    static int[][] erodeGrid(double amount){
        int upperbound = (int)Math.ceil(amount);
        double amountsq = (amount + 0.5) * (amount + 0.5);// adding 0.5 due to the erosion being from the center of the pixel
        //(re)-initialize the edited_grid
        int[][] edited_grid = new int[grid.length][grid[0].length];
        for(int x = 0; x < grid.length; x++){
            for(int y = 0; y < grid[0].length; y++){
                edited_grid[x][y] = grid[x][y];
            }
        }
        if(amount >= 0.2){// don't bother with the erosion if the amount is too small to change anything
            //iterate through every point
            for(int x = 0; x < grid.length; x++){
                for(int y = 0; y < grid[0].length; y++){
                    if(grid[x][y] >= 1){
                        //if the value of the grid point can be changed, iterate through its neighboring points
                        int[][] subdivisions = new int[5][5];
                        for(int xx = - upperbound; xx <= upperbound; xx++){
                            for(int yy = - upperbound; yy <= upperbound; yy++){
                                if(getGrid(x+xx, y+yy) == 0){
                                    /*if the point at (x+xx,y+yy) is different from 
                                    the one at (x,y), subdivide the point at (x,y)
                                    into a 5x5 grid and see how many of these grid
                                    points are within the erosion or dilation distance
                                    */
                                    for(int xxx = 0; xxx < 5; xxx ++){
                                        for(int yyy = 0; yyy < 5; yyy ++){
                                            double x_subdist = ((double)xxx)*0.2 - 0.4;
                                            double y_subdist = ((double)yyy)*0.2 - 0.4;
                                            double xdist = ((double)xx) + x_subdist;
                                            double ydist = ((double)yy) + y_subdist;
                                            double distancesq = xdist*xdist + ydist*ydist;
                                            if(distancesq < amountsq){
                                                subdivisions[xxx][yyy] = 1;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        /*compute the number of subdivisions at (x,y) that should change.
                        If they are in the majority, flip the pixel.
                        */
                        int sum = 0;
                        for(int i = 0; i < 5; i++){
                            for(int ii = 0; ii < 5; ii++){
                                sum += subdivisions[i][ii];
                            }
                        }
                        if(sum > 12){
                            edited_grid[x][y] = 0;
                        }
                        
                    }
                }
            }
        }
        return edited_grid;
    }
    
    /**
     * Counts the number of pixels that are switched on in a given grid
     * @param de_grid A grid where '1' indicates the pixel is on and '0' indicates the pixel is off.
     * @return 
     */
    static int getSum(int[][] de_grid){
        int runningcount = 0;
        for(int x = 0; x < de_grid.length; x++){
            for(int y = 0; y < de_grid[0].length; y++){
                runningcount += de_grid[x][y];
            }
        }
        return runningcount;
    }
    
    /**
     * Gets the value of the {@link #grid} at (x,y). If (x,y) is outside of the grid, 0 will be returned instead.
     * @param x The x-index of interest.
     * @param y The y-index of interest.
     * @return The value of the grid[x][y] if x and y are in bounds. Zero otherwise.
     */
    static int getGrid(int x, int y){
        if(x > 0 && y > 0 && x < grid.length && y < grid[0].length){
            return grid[x][y];
        }
        else{
            return 0;
        }
    }
    
    /*
    TESTING AND DEBUGGING PURPOSES ONLY
    */
    static void exportEditedGrid(String filename, int[][] edited_grid) throws IOException{
        File f = new File(filename);
        f.createNewFile();
        PrintWriter pw = new PrintWriter(f);
        pw.print("#x,y,switch");
        for(int x = 0; x < edited_grid.length; x++){
            for(int y = 0; y < edited_grid[0].length; y++){
                pw.print("\n" + x + "," + y + "," + edited_grid[x][y]);
            }
        }
        pw.close();
    }
}
