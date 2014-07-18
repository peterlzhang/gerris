/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001-2012 National Institute of Water and Atmospheric Research
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.  
 */
/*! \file
 * \brief Boundary conditions.
 */

#include <stdlib.h>
#include <math.h>

#include "config.h"
#include "boundary.h"
#include "simulation.h"
#include "adaptive.h"
#include "vof.h"

static FttVector rpos[FTT_NEIGHBORS] = {
#if FTT_2D
  {1.,0.,0.}, {-1.,0.,0.}, {0.,1.,0.}, {0.,-1.,0.}
#else  /* FTT_3D */
  {1.,0.,0.}, {-1.,0.,0.}, {0.,1.,0.}, {0.,-1.,0.}, {0.,0.,1.}, {0.,0.,-1.}
#endif /* FTT_3D */
};

/**
 * Default symmetry boundary condition.
 * \beginobject{GfsBc}
 */

static void symmetry (FttCellFace * f, GfsBc * b)
{
  if (b->v->component == f->d/2 && !b->v->even)
    GFS_VALUE (f->cell, b->v) = - GFS_VALUE (f->neighbor, b->v);
  else
    GFS_VALUE (f->cell, b->v) =   GFS_VALUE (f->neighbor, b->v);
}

static void set_stencil_neighbor (FttCellFace * f, GfsBc * b, gdouble w)
{
  GFS_DOUBLE_TO_POINTER (GFS_VALUE (f->cell, b->lp->neighbor)) = f->neighbor;
  GFS_VALUE (f->cell, b->lp->neighborw) = w;
}

static void symmetry_stencil (FttCellFace * f, GfsBc * b)
{
  set_stencil_neighbor (f, b, (b->v->component == f->d/2 && !b->v->even) ? -1. : 1.);
}

static void face_symmetry (FttCellFace * f, GfsBc * b)
{
  if (b->v->component == f->d/2 && !b->v->even)
    GFS_STATE (f->cell)->f[f->d].v = 
      GFS_STATE (f->neighbor)->f[FTT_OPPOSITE_DIRECTION (f->d)].v = 0.;
  else if (GFS_IS_VARIABLE_TRACER_VOF (b->v))
    GFS_STATE (f->cell)->f[f->d].v = GFS_VALUE (f->neighbor, b->v);
  else
    GFS_STATE (f->cell)->f[f->d].v = 
      GFS_STATE (f->neighbor)->f[FTT_OPPOSITE_DIRECTION (f->d)].v;
}

static void bc_write (GtsObject * o, FILE * fp)
{
  g_assert (GFS_BC (o)->v);
  fprintf (fp, "%s %s", o->klass->info.name, GFS_BC (o)->v->name);
}

static void bc_read (GtsObject ** o, GtsFile * fp)
{
  GfsBc * bc = GFS_BC (*o);

  g_assert (bc->b);

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (klass)");
    return;
  }
  gts_file_next_token (fp);

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (v)");
    return;
  }
  bc->v = gfs_variable_from_name (gfs_box_domain (bc->b->box)->variables, 
				  fp->token->str);
  if (bc->v == NULL)
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
  else
    gts_file_next_token (fp);
}

static void gfs_bc_class_init (GtsObjectClass * klass)
{
  klass->write = bc_write;
  klass->read =  bc_read;
}

static void gfs_bc_init (GfsBc * object)
{
  object->bc =                     (FttFaceTraverseFunc) symmetry;
  object->homogeneous_bc =         (FttFaceTraverseFunc) symmetry;
  object->homogeneous_bc_stencil = (FttFaceTraverseFunc) symmetry_stencil;
  object->face_bc =                (FttFaceTraverseFunc) face_symmetry;
}

GfsBcClass * gfs_bc_class (void)
{
  static GfsBcClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_bc_info = {
      "GfsBc",
      sizeof (GfsBc),
      sizeof (GfsBcClass),
      (GtsObjectClassInitFunc) gfs_bc_class_init,
      (GtsObjectInitFunc) gfs_bc_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_object_class ()),
				  &gfs_bc_info);
  }

  return klass;
}

/**
 * gfs_bc_new:
 * @k: a #GfsBcClass.
 * @v: a variable associated with the BC or %NULL.
 * @extra:
 *
 * Returns: a new #GfsBc.
 */
GfsBc * gfs_bc_new (GfsBcClass * k, GfsVariable * v, gboolean extra)
{
  GfsBc * b;

  g_return_val_if_fail (k != NULL, NULL);

  b = GFS_BC (gts_object_new (GTS_OBJECT_CLASS (k)));
  if (v)
    gfs_object_simulation_set (b, v->domain);
  b->v = v;
  b->extra = extra;

  return b;
}

/** \endobject{GfsBc} */

/**
 * Generic class for a boundary condition taking a parameter.
 * \beginobject{GfsBcValue}
 */

static void bc_value_write (GtsObject * o, FILE * fp)
{
  if (GTS_OBJECT_CLASS (gfs_bc_value_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_bc_value_class ())->parent_class->write) 
      (o, fp);
  if (GFS_BC_VALUE (o)->val)
    gfs_function_write (GFS_BC_VALUE (o)->val, fp);
}

static void bc_value_read (GtsObject ** o, GtsFile * fp)
{
  GfsBcValue * bc = GFS_BC_VALUE (*o);

  if (GTS_OBJECT_CLASS (gfs_bc_value_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_bc_value_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  gfs_function_read (GFS_BC_VALUE (*o)->val, gfs_box_domain (GFS_BC (bc)->b->box), fp);
}

static void bc_value_destroy (GtsObject * o)
{
  gts_object_destroy (GTS_OBJECT (GFS_BC_VALUE (o)->val));
  (* GTS_OBJECT_CLASS (gfs_bc_value_class ())->parent_class->destroy) (o);
}

static void gfs_bc_value_class_init (GtsObjectClass * klass)
{
  klass->write   = bc_value_write;
  klass->read    = bc_value_read;
  klass->destroy = bc_value_destroy;
}

static void gfs_bc_value_init (GfsBcValue * bc)
{
  bc->val = gfs_function_new (gfs_function_class (), 0.);
}

GfsBcClass * gfs_bc_value_class (void)
{
  static GfsBcClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_bc_value_info = {
      "GfsBcValue",
      sizeof (GfsBcValue),
      sizeof (GfsBcClass),
      (GtsObjectClassInitFunc) gfs_bc_value_class_init,
      (GtsObjectInitFunc) gfs_bc_value_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_bc_class ()),
				  &gfs_bc_value_info);
  }

  return klass;
}

static GfsBc * gfs_bc_value_new (GfsBcClass * k,
				 GfsVariable * v,
				 GfsFunction * val,
				 gboolean extra)
{
  GfsBcValue * bc = GFS_BC_VALUE (gfs_bc_new (k, v, extra));

  if (val != NULL) {
    gts_object_destroy (GTS_OBJECT (bc->val));
    bc->val = val;
  }

  return GFS_BC (bc);
}

/** \endobject{GfsBcValue} */

/**
 * Dirichlet boundary condition.
 * \beginobject{GfsBcDirichlet}
 */

static void dirichlet (FttCellFace * f, GfsBc * b) // FttCellFace structure stores cell, neighbor, d.
{ static int count = 0;
  count += 1;
  GFS_VALUE (f->cell, b->v) = 
    2.*gfs_function_face_value (GFS_BC_VALUE (b)->val, f)
    - GFS_VALUE (f->neighbor, b->v);
  //printf("Dirichlet count = %d, val at interior point = %f\n",count,GFS_VALUE (f->neighbor, b->v));
}

static void dirichlet_vof (FttCellFace * f, GfsBc * b)
{
  GFS_VALUE (f->cell, b->v) = gfs_function_face_value (GFS_BC_VALUE (b)->val, f);
}

static void homogeneous_dirichlet (FttCellFace * f, GfsBc * b)
{
  GFS_VALUE (f->cell, b->v) = - GFS_VALUE (f->neighbor, b->v);
}

static void homogeneous_dirichlet_stencil (FttCellFace * f, GfsBc * b)
{
  set_stencil_neighbor (f, b, -1.);
}

static void face_dirichlet (FttCellFace * f, GfsBc * b)
{
  GFS_STATE (f->cell)->f[f->d].v = GFS_STATE (f->neighbor)->f[FTT_OPPOSITE_DIRECTION (f->d)].v = 
    gfs_function_face_value (GFS_BC_VALUE (b)->val, f);
}

static void bc_dirichlet_read (GtsObject ** o, GtsFile * fp)
{
  GfsBc * bc = GFS_BC (*o);

  if (GTS_OBJECT_CLASS (gfs_bc_dirichlet_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_bc_dirichlet_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  gfs_function_set_units (GFS_BC_VALUE (bc)->val, bc->v->units);
  if (GFS_IS_VARIABLE_TRACER_VOF (bc->v))
    bc->bc = (FttFaceTraverseFunc) dirichlet_vof;
}

static void gfs_bc_dirichlet_init (GfsBc * object)
{
  object->bc =                     (FttFaceTraverseFunc) dirichlet;
  object->homogeneous_bc =         (FttFaceTraverseFunc) homogeneous_dirichlet;
  object->homogeneous_bc_stencil = (FttFaceTraverseFunc) homogeneous_dirichlet_stencil;
  object->face_bc =                (FttFaceTraverseFunc) face_dirichlet;
}

static void gfs_bc_dirichlet_class_init (GtsObjectClass * klass)
{
  klass->read = bc_dirichlet_read;
}

GfsBcClass * gfs_bc_dirichlet_class (void)
{
  static GfsBcClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_bc_dirichlet_info = {
      "GfsBcDirichlet",
      sizeof (GfsBcValue),
      sizeof (GfsBcClass),
      (GtsObjectClassInitFunc) gfs_bc_dirichlet_class_init,
      (GtsObjectInitFunc) gfs_bc_dirichlet_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_bc_value_class ()),
				  &gfs_bc_dirichlet_info);
  }

  return klass;
}

/** \endobject{GfsBcDirichlet} */

/**
 * Neumann boundary condition.
 * \beginobject{GfsBcNeumann}
 */

static void neumann (FttCellFace * f, GfsBc * b)
{
  GFS_VALUE (f->cell, b->v) = 
    GFS_VALUE (f->neighbor, b->v) +
    gfs_function_face_value (GFS_BC_VALUE (b)->val, f)
    *ftt_cell_size (f->cell);
}

static void homogeneous_neumann (FttCellFace * f, GfsBc * b)
{
  GFS_VALUE (f->cell, b->v) = GFS_VALUE (f->neighbor, b->v);
}

static void homogeneous_neumann_stencil (FttCellFace * f, GfsBc * b)
{
  set_stencil_neighbor (f, b, 1.);
}

static void face_neumann (FttCellFace * f, GfsBc * b)
{
  GFS_STATE (f->cell)->f[f->d].v = 
    GFS_VALUE (f->neighbor, b->v) +
    gfs_function_face_value (GFS_BC_VALUE (b)->val, f)
    *ftt_cell_size (f->cell)/2.;
}

static void bc_neumann_read (GtsObject ** o, GtsFile * fp)
{
  GfsBc * bc = GFS_BC (*o);

  if (GTS_OBJECT_CLASS (gfs_bc_neumann_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_bc_neumann_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  gfs_function_set_units (GFS_BC_VALUE (bc)->val, bc->v->units - 1.);
}

static void gfs_bc_neumann_init (GfsBc * object)
{
  object->bc =                     (FttFaceTraverseFunc) neumann;
  object->homogeneous_bc =         (FttFaceTraverseFunc) homogeneous_neumann;
  object->homogeneous_bc_stencil = (FttFaceTraverseFunc) homogeneous_neumann_stencil;
  object->face_bc =                (FttFaceTraverseFunc) face_neumann;
}

static void gfs_bc_neumann_class_init (GtsObjectClass * klass)
{
  klass->read = bc_neumann_read;
}

GfsBcClass * gfs_bc_neumann_class (void)
{
  static GfsBcClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_bc_neumann_info = {
      "GfsBcNeumann",
      sizeof (GfsBcValue),
      sizeof (GfsBcClass),
      (GtsObjectClassInitFunc) gfs_bc_neumann_class_init,
      (GtsObjectInitFunc) gfs_bc_neumann_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_bc_value_class ()),
				  &gfs_bc_neumann_info);
  }

  return klass;
}

/** \endobject{GfsBcNeumann} */

/**
 * Contact angle boundary condition.
 * \beginobject{GfsBcAngle}
 */

static void bc_angle_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_bc_angle_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (!GFS_IS_VARIABLE_TRACER_VOF_HEIGHT (GFS_BC (*o)->v))
    gts_file_error (fp, "expecting a GfsVariableTracerVOFHeight");
  gfs_function_set_units (GFS_BC_VALUE (*o)->val, 0.);
}

static void gfs_bc_angle_init (GfsBc * object)
{
  /* use zero for Neumann condition for the VOF tracer */
  object->bc = (FttFaceTraverseFunc) homogeneous_neumann;
}

static void gfs_bc_angle_class_init (GtsObjectClass * klass)
{
  klass->read = bc_angle_read;
}

GfsBcClass * gfs_bc_angle_class (void)
{
  static GfsBcClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsBcAngle",
      sizeof (GfsBcValue),
      sizeof (GfsBcClass),
      (GtsObjectClassInitFunc) gfs_bc_angle_class_init,
      (GtsObjectInitFunc) gfs_bc_angle_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_bc_neumann_class ()), &info);
  }

  return klass;
}

/** \endobject{GfsBcAngle} */

/*----------------------------------------------------------------------------------------------------------*/
/**
 * Joseph boundary condition.
 * \beginobject{GfsBcJoseph}
 */

/* compute tangential gradient of varaible @v along the boundary
 * 2nd order accurate velocity gradient
 * takes a one sided difference cell on one side doesn't exist
 */

static gdouble boundary_tangential_gradient(FttCellFace * f, GfsVariable * v)
{
  gdouble dudt;
  FttComponent c = f->d % 2;
  FttComponent oc = FTT_ORTHOGONAL_COMPONENT(c);
  FttDirection nd = 2*oc;
  FttCell * nlg = ftt_cell_neighbor (f->cell, nd), * nrg = ftt_cell_neighbor (f->cell, nd + 1);
  FttCell * nl = ftt_cell_neighbor (f->neighbor,nd), * nr = ftt_cell_neighbor (f->neighbor, nd + 1);

  gdouble h = ftt_cell_size(f->cell);
  FttVector pos, posl, posr;
  
  if (!nlg) { // at corner, no left ghost cell
    FttCell * nrrg = ftt_cell_neighbor (nrg, nd+1);
    FttCell * nrr = ftt_cell_neighbor (nr, nd+1);
    
    gdouble vb = (GFS_VALUE(f->cell,v)+GFS_VALUE(f->neighbor,v))/2.;
    gdouble vr = (GFS_VALUE(nrg,v)+GFS_VALUE(nr,v))/2.;
    gdouble vrr = (GFS_VALUE(nrrg,v)+GFS_VALUE(nrr,v))/2.;

    dudt = (-3.*vb+4.*vr-vrr)/(2.*h);
  }
  else if (!nrg) { // at corner, no left ghost cell
    FttCell * nllg = ftt_cell_neighbor (nlg, nd);
    FttCell * nll = ftt_cell_neighbor (nl, nd);
    
    gdouble vb = (GFS_VALUE(f->cell,v)+GFS_VALUE(f->neighbor,v))/2.;
    gdouble vl = (GFS_VALUE(nlg,v)+GFS_VALUE(nl,v))/2.;
    gdouble vll = (GFS_VALUE(nllg,v)+GFS_VALUE(nll,v))/2.;

    dudt = (3.*vb-4.*vl+vll)/(2.*h);
  }
  else { // left and right neighbors exist
    gdouble vl = (GFS_VALUE(nlg,v)+GFS_VALUE(nl,v))/2.;
    gdouble vr = (GFS_VALUE(nrg,v)+GFS_VALUE(nr,v))/2.;
    
    dudt = (vr-vl)/(2.*h);
  }

  return dudt;
}

static gboolean is_interfacial (FttCell * cell, gpointer data)
{
  GfsVariable * f = data;
  return (GFS_VALUE (cell, f) > 0. && GFS_VALUE (cell, f) < 1.);
}

/* find_tcp locates the triple contact point at the face defined 
 * between f->cell and f->neighbor.
 * vof tracer is defined by the variable b->vofv
 */
static FttVector find_tcp(FttCellFace * f, GfsBc *b)
{
  g_assert (is_interfacial(f->neighbor,b->vofv));

  guint i;
  gdouble h = ftt_cell_size(f->cell);
  FttVector p1, p2, w1, w2, tcp, interfacenorm, pos, posn;
  GfsVariableTracerVOF * t = GFS_VARIABLE_TRACER_VOF(b->vofv);
  FttVector q[FTT_DIMENSION*(FTT_DIMENSION - 1) + 1];
  guint ndim = gfs_vof_facet (f->neighbor, t, q, &interfacenorm);
  for (i = 0; i < ndim-1; i++ ) {
    p1.x = q[i].x; p1.y = q[i].y; p1.z = q[i].z;
    p2.x = q[i + 1].x; p2.y = q[i + 1].y; p2.z = q[i + 1].z;
  }

  /* Wall end points */
  ftt_cell_pos(f->cell,&pos);
  ftt_cell_pos(f->neighbor,&posn);
  if (f->d/2 == 0) {
    w1.x = (pos.x+posn.x)/2.; w1.y = pos.y + h/2.;
    w2.x = (pos.x+posn.x)/2.; w2.y = pos.y - h/2.;
  }
  else if (f->d/2 == 1) {
    w1.y = (pos.y+posn.y)/2.; w1.x = pos.x + h/2.;
    w2.y = (pos.y+posn.y)/2.; w2.x = pos.x - h/2.;    
  }

//  printf("cell(%f,%f), h = %f, w1(%f,%f), w2(%f,%f)\n",pos.x,pos.y,h,w1.x,w1.y,w2.x,w2.y);

  tcp.x = ((w1.x*w2.y-w1.y*w2.x)*(p1.x-p2.x)-(w1.x-w2.x)*(p1.x*p2.y-p1.y*p2.x))/
              ((w1.x-w2.x)*(p1.y-p2.y)-(w1.y-w2.y)*(p1.x-p2.x));
  tcp.y = ((w1.x*w2.y-w1.y*w2.x)*(p1.y-p2.y)-(w1.y-w2.y)*(p1.x*p2.y-p1.y*p2.x))/
              ((w1.x-w2.x)*(p1.y-p2.y)-(w1.y-w2.y)*(p1.x-p2.x));

  return tcp;
}

typedef struct {
  gdouble phi,mua,mub,Lsa,Lsb,alpha,U,r,theta;
} stokes_params;

static gdouble stokes_stress_integral(gdouble (*f)(stokes_params,gint),gdouble a,gdouble b,gint n,stokes_params data,gint option)
{
  gdouble h,sum,x0,term;
  gint k;
  stokes_params sp1, sp2, sp3;

  if (option == 0 || option == 1) {
    data.theta = M_PI;
  }
  else if (option == 2 || option == 3) {
    data.theta = 0;
  }

  sp1 = data;
  sp2 = data;
  sp3 = data;

  h = (b - a)/n;
  sum = 0.0;

  for(k = 0;k < n;k++)
  {
    x0 = a + k*h;
    sp1.r = x0;
    sp2.r = x0+h/2.;
    sp3.r = x0+h;
    term = h*(f(sp1,option)+4*f(sp2,option)+f(sp3,option))/6;
    sum += term;
  }

  return sum;
}

static gdouble stokes_v_gradient(stokes_params data,gint option)
{
  gdouble phi = data.phi;
  gdouble mua = data.mua;
  gdouble mub = data.mub;
  gdouble Lsa = data.Lsa;
  gdouble Lsb = data.Lsb;
  gdouble alpha = data.alpha;
  gdouble U = data.U;
  gdouble r = data.r;
  gdouble theta = data.theta;

  gdouble c = cos(phi),s = sin(phi), pi=M_PI;
  gdouble p1, p2, p3, p4;

  p1 = s*s-phi*phi;
  p2 = phi-c*s;
  p3 = (phi-pi*c*c);
  p4 = phi-pi;

  gdouble D1, D2, D3, D4, D5, D6, D7, D8, D9, N1, N2, N3, N4, N5, N6, aa, ba, ca, da, ab, bb, cb, db,
    aaD1, aaD2, aaD3, aaD4, aaD5, aaD6, aaD7, aaD8, aaD9,
    baD1, baD2, baD3, baD4, baD5, baD6, baD7, baD8, baD9,
    caD1, caD2, caD3, caD4, caD5, caD6, caD7, caD8, caD9,
    daD1, daD2, daD3, daD4, daD5, daD6, daD7, daD8, daD9,
    abD1, abD2, abD3, abD4, abD5, abD6, abD7, abD8, abD9,
    cbD1, cbD2, cbD3, cbD4, cbD5, cbD6, cbD7, cbD8, cbD9,
    dbD1, dbD2, dbD3, dbD4, dbD5, dbD6, dbD7, dbD8, dbD9,
    aa_numer, aa_denom, ba_numer, ba_denom, ca_numer, ca_denom, da_numer, da_denom,
    ab_numer, ab_denom, cb_numer, cb_denom, db_numer, db_denom,
    daadr, daardr, dbadr, dbardr, dcadr, dcardr, ddadr, ddardr,
    dabdr, dabrdr, dbbdr, dbbrdr, dcbdr, dcbrdr, ddbdr, ddbrdr,
    psia,psib,dpsidthetaa,dpsidthetab,dpsidra,dpsidrb,
    dpsidtheta2a,dpsidthetadra,dpsidtheta2b,dpsidthetadrb,dvrdthetaa,dvrdra,dvrdthetab,dvrdrb;

  D1 = alpha*(phi-c*s)*(pi-phi+c*s);
  D2 = (s*s-phi*phi)*(pi-phi+c*s);
  D3 = (phi-c*s)*(2.*phi*pi-phi*phi-pi*pi+s*s);
  D4 = alpha*s*s*(Lsa*phi-Lsb*phi+Lsb*pi-Lsa*c*s+Lsb*c*s);
  D5 = Lsa*s*s*(s*s-phi*phi)+Lsb*(phi*phi-2.*phi*c*s-pi*phi+c*c*s*s+pi*c*s);
  D6 = Lsa*(phi*phi-2.*phi*c*s-pi*phi+c*c*s*s+pi*c*s)+Lsb*s*s*(2.*phi*pi-phi*phi-pi*pi+s*s);
  D7 = Lsa*Lsb*alpha*s*s*s*s;
  D8 = Lsa*Lsb*s*s*(phi-c*s);
  D9 = Lsa*Lsb*s*s*(pi-phi+c*s);

    // aa ----------------------------------------------------------
  N1 = -16*U*alpha*p3*p2;
  N2 = -8*U*p3*p1;
  N3 = -8*U*(p4*(phi*phi-phi*pi*c*c+pi*c*s)-s*s*p3);
  N4 = -32*Lsb*U*alpha*s*s*(pi*s*s+phi-pi);
  N5 = 16*Lsb*U*p3*p2;
  N6 = 16*Lsa*U*p4*(phi*c*c-c*s)+16*Lsb*U*s*s*(phi*pi-phi*phi+pi*c*s);

  aaD1 = -16*D1;
  aaD2 = -8*D2;
  aaD3 = -8*D3;
  aaD4 = -32*D4;
  aaD5 = -16*D5;
  aaD6 = -16*D6;
  aaD7 = -64*D7;
  aaD8 = 32*D8;
  aaD9 = 32*D9;

  aa_numer = mua*mub*r*r*N1+mua*r*r*N2+mub*r*r*N3+mua*mub*r*N4+mua*r*N5+mub*r*N6;
  aa_denom = mua*mub*(r*r*aaD1+r*aaD4+aaD7)+mua*(r*r*aaD2+r*aaD5+aaD8)+
      mub*(r*r*aaD3+r*aaD6+aaD9);

  aa = aa_numer/aa_denom;

  daadr = (aa_denom*(2*r*mua*mub*N1+2*r*mua*N2+2*r*mub*N3+mua*mub*N4+mua*N5+mub*N6)-
      aa_numer*(2*r*mua*mub*aaD1+2*r*mua*aaD2+2*r*mub*aaD3+mua*mub*aaD4+mua*aaD5+mub*aaD6))/
      (aa_denom*aa_denom);

  daardr = r*daadr+aa;
    // ba ----------------------------------------------------------
  N1 = 2*U*alpha*pi*c*s*p2;
  N2 = U*pi*c*s*p1;
  N3 = U*pi*s*(phi*phi*c-c*s*s+pi*s-phi*pi*c);
  N4 = 4*Lsb*U*alpha*pi*c*s*s*s;
  N5 = -2*Lsb*U*pi*c*s*p2;
  N6 = 2*U*pi*(Lsb*s*s*s*s+Lsa*s*s-Lsa*phi*c*s);

  baD1 = 2*D1;
  baD2 = 1*D2;
  baD3 = 1*D3;
  baD4 = 4*D4;
  baD5 = 2*D5;
  baD6 = 2*D6;
  baD7 = 8*D7;
  baD8 = -4*D8;
  baD9 = -4*D9;

  ba_numer = mua*mub*r*r*N1+mua*r*r*N2+mub*r*r*N3+mua*mub*r*N4+mua*r*N5+mub*r*N6;
  ba_denom = mua*mub*(r*r*baD1+r*baD4+baD7)+mua*(r*r*baD2+r*baD5+baD8)+
      mub*(r*r*baD3+r*baD6+baD9);

  ba = ba_numer/ba_denom;

  dbadr = (ba_denom*(2*r*mua*mub*N1+2*r*mua*N2+2*r*mub*N3+mua*mub*N4+mua*N5+mub*N6)-
      ba_numer*(2*r*mua*mub*baD1+2*r*mua*baD2+2*r*mub*baD3+mua*mub*baD4+mua*baD5+mub*baD6))/
      (ba_denom*ba_denom);

  dbardr = r*dbadr+ba;
    // ca ----------------------------------------------------------
  N1 = -8*U*alpha*s*s*p2;
  N2 = -4*U*s*s*p1;
  N3 = 4*U*s*s*(pi*phi-phi*phi+s*s);
  N4 = -16*Lsb*U*alpha*s*s*s*s;
  N5 = 8*Lsb*U*s*s*p2;
  N6 = 8*Lsb*U*s*s*(c*s-phi+pi);

  caD1 = 8*D1;
  caD2 = 4*D2;
  caD3 = 4*D3;
  caD4 = 16*D4;
  caD5 = 8*D5;
  caD6 = 8*D6;
  caD7 = 32*D7;
  caD8 = -16*D8;
  caD9 = -16*D9;

  ca_numer = mua*mub*r*r*N1+mua*r*r*N2+mub*r*r*N3+mua*mub*r*N4+mua*r*N5+mub*r*N6;
  ca_denom = mua*mub*(r*r*caD1+r*caD4+caD7)+mua*(r*r*caD2+r*caD5+caD8)+
      mub*(r*r*caD3+r*caD6+caD9);

  ca = ca_numer/ca_denom;

  dcadr = (ca_denom*(2*r*mua*mub*N1+2*r*mua*N2+2*r*mub*N3+mua*mub*N4+mua*N5+mub*N6)-
      ca_numer*(2*r*mua*mub*caD1+2*r*mua*caD2+2*r*mub*caD3+mua*mub*caD4+mua*caD5+mub*caD6))/
      (ca_denom*ca_denom);

  dcardr = r*dcadr+ca;
    // da ----------------------------------------------------------
  N1 = -2*U*alpha*c*s*p2;
  N2 = -U*c*s*p1;
  N3 = U*c*s*p1-U*s*pi*(s-phi*c);
  N4 = -4*Lsb*U*alpha*c*s*s*s;
  N5 = 2*Lsb*U*c*s*p2;
  N6 = -2*U*(Lsb*s*s*s*s+Lsa*s*s-Lsa*phi*c*s);

  daD1 = 2*D1;
  daD2 = 1*D2;
  daD3 = 1*D3;
  daD4 = 4*D4;
  daD5 = 2*D5;
  daD6 = 2*D6;
  daD7 = 8*D7;
  daD8 = -4*D8;
  daD9 = -4*D9;

  da_numer = mua*mub*r*r*N1+mua*r*r*N2+mub*r*r*N3+mua*mub*r*N4+mua*r*N5+mub*r*N6;
  da_denom = mua*mub*(r*r*daD1+r*daD4+daD7)+mua*(r*r*daD2+r*daD5+daD8)+
      mub*(r*r*daD3+r*daD6+daD9);

  da = da_numer/da_denom;

  ddadr = (da_denom*(2*r*mua*mub*N1+2*r*mua*N2+2*r*mub*N3+mua*mub*N4+mua*N5+mub*N6)-
      da_numer*(2*r*mua*mub*daD1+2*r*mua*daD2+2*r*mub*daD3+mua*mub*daD4+mua*daD5+mub*daD6))/
      (da_denom*da_denom);

  ddardr = r*ddadr+da;
    // ab ----------------------------------------------------------
  N1 = -8*U*alpha*phi*(pi-phi+c*s);
  N2 = 4*U*phi*(phi*pi+s*s-phi*phi+pi*c*s);
  N3 = -4*U*phi*(s*s-phi*phi+2*phi*pi-pi*pi);
  N4 = -16*Lsa*U*alpha*phi*s*s;
  N5 = 8*Lsa*U*phi*phi*s*s+8*Lsb*U*phi*(pi*c*c-phi*c*c+c*s);
  N6 = 8*Lsa*U*phi*(pi-phi+c*s);

  abD1 = 8*D1;
  abD2 = 4*D2;
  abD3 = 4*D3;
  abD4 = 16*D4;
  abD5 = 8*D5;
  abD6 = 8*D6;
  abD7 = 32*D7;
  abD8 = -16*D8;
  abD9 = -16*D9;

  ab_numer = mua*mub*r*r*N1+mua*r*r*N2+mub*r*r*N3+mua*mub*r*N4+mua*r*N5+mub*r*N6;
  ab_denom = mua*mub*(r*r*abD1+r*abD4+abD7)+mua*(r*r*abD2+r*abD5+abD8)+
      mub*(r*r*abD3+r*abD6+abD9);

  ab = ab_numer/ab_denom;

  dabdr = (ab_denom*(2*r*mua*mub*N1+2*r*mua*N2+2*r*mub*N3+mua*mub*N4+mua*N5+mub*N6)-
      ab_numer*(2*r*mua*mub*abD1+2*r*mua*abD2+2*r*mub*abD3+mua*mub*abD4+mua*abD5+mub*abD6))/
      (ab_denom*ab_denom);

  dabrdr = r*dabdr+ab;
    // bb ----------------------------------------------------------

  bb = 0;

  dbbdr = 0;
  dbbrdr = 0;
    // cb ----------------------------------------------------------
  N1 = 8*U*alpha*s*s*(c*s-phi+pi);
  N2 = -4*U*s*s*(s*s+pi*phi-phi*phi);
  N3 = 4*U*s*s*(2*phi*pi-phi*phi-pi*pi+s*s);
  N4 = 16*Lsa*U*alpha*s*s*s*s;
  N5 = -8*Lsa*U*s*s*p2;
  N6 = -8*Lsa*U*s*s*(c*s-phi+pi);

  cbD1 = 8*D1;
  cbD2 = 4*D2;
  cbD3 = 4*D3;
  cbD4 = 16*D4;
  cbD5 = 8*D5;
  cbD6 = 8*D6;
  cbD7 = 32*D7;
  cbD8 = -16*D8;
  cbD9 = -16*D9;

  cb_numer = mua*mub*r*r*N1+mua*r*r*N2+mub*r*r*N3+mua*mub*r*N4+mua*r*N5+mub*r*N6;
  cb_denom = mua*mub*(r*r*cbD1+r*cbD4+cbD7)+mua*(r*r*cbD2+r*cbD5+cbD8)+
      mub*(r*r*cbD3+r*cbD6+cbD9);

  cb = cb_numer/cb_denom;

  dcbdr = (cb_denom*(2*r*mua*mub*N1+2*r*mua*N2+2*r*mub*N3+mua*mub*N4+mua*N5+mub*N6)-
      cb_numer*(2*r*mua*mub*cbD1+2*r*mua*cbD2+2*r*mub*cbD3+mua*mub*cbD4+mua*cbD5+mub*cbD6))/
      (cb_denom*cb_denom);

  dcbrdr = r*dcbdr+cb;
    // db ----------------------------------------------------------
  N1 = 2*U*alpha*c*s*(pi-phi+c*s);
  N2 = U*s*(phi*phi*c-c*s*s-pi*s-phi*pi*c);
  N3 = U*c*s*(2*phi*pi-phi*phi-pi*pi+s*s);
  N4 = 4*Lsa*U*alpha*c*s*s*s;
  N5 = -2*U*(Lsa*s*s*s*s+Lsb*s*s-Lsb*phi*c*s+Lsb*pi*c*s);
  N6 = -2*Lsa*U*c*s*(pi-phi+c*s);

  dbD1 = 2*D1;
  dbD2 = 1*D2;
  dbD3 = 1*D3;
  dbD4 = 4*D4;
  dbD5 = 2*D5;
  dbD6 = 2*D6;
  dbD7 = 8*D7;
  dbD8 = -4*D8;
  dbD9 = -4*D9;

  db_numer = mua*mub*r*r*N1+mua*r*r*N2+mub*r*r*N3+mua*mub*r*N4+mua*r*N5+mub*r*N6;
  db_denom = mua*mub*(r*r*dbD1+r*dbD4+dbD7)+mua*(r*r*dbD2+r*dbD5+dbD8)+
      mub*(r*r*dbD3+r*dbD6+dbD9);

  db = db_numer/db_denom;

  ddbdr = (db_denom*(2*r*mua*mub*N1+2*r*mua*N2+2*r*mub*N3+mua*mub*N4+mua*N5+mub*N6)-
      db_numer*(2*r*mua*mub*dbD1+2*r*mua*dbD2+2*r*mub*dbD3+mua*mub*dbD4+mua*dbD5+mub*dbD6))/
      (db_denom*db_denom);

  ddbrdr = r*ddbdr+db;
  /////////////////////////////////
  psia = r*(aa*sin(theta)+ba*cos(theta)+ca*theta*sin(theta)+da*theta*cos(theta));
  psib = r*(ab*sin(theta)+bb*cos(theta)+cb*theta*sin(theta)+db*theta*cos(theta));

  dpsidthetaa = r*(aa*cos(theta)-ba*sin(theta)+ca*sin(theta)+
      ca*theta*cos(theta)-da*theta*sin(theta)+da*cos(theta));
  dpsidthetab = r*(ab*cos(theta)-bb*sin(theta)+cb*sin(theta)+
      cb*theta*cos(theta)-db*theta*sin(theta)+db*cos(theta));
  dpsidra = daardr*sin(theta)+dbardr*cos(theta)+dcardr*theta*sin(theta)+ddardr*theta*cos(theta);
  dpsidrb = dabrdr*sin(theta)+dbbrdr*cos(theta)+dcbrdr*theta*sin(theta)+ddbrdr*theta*cos(theta);


  dpsidtheta2a = r*(-aa*sin(theta)-ba*cos(theta)+2*ca*cos(theta)-
      ca*theta*sin(theta)-2*da*sin(theta)-da*theta*cos(theta));
  dpsidthetadra = daardr*cos(theta)-dbardr*sin(theta)+
      dcardr*sin(theta)+dcardr*theta*cos(theta)-ddardr*theta*sin(theta)+ddardr*cos(theta);

  dpsidtheta2b = r*(-ab*sin(theta)-bb*cos(theta)+2*cb*cos(theta)-
      cb*theta*sin(theta)-2*db*sin(theta)-db*theta*cos(theta));
  dpsidthetadrb = dabrdr*cos(theta)-dbbrdr*sin(theta)+
      dcbrdr*sin(theta)+dcbrdr*theta*cos(theta)-ddbrdr*theta*sin(theta)+ddbrdr*cos(theta);


  dvrdthetaa = -(1/(r*r))*dpsidtheta2a;
  dvrdra = (1/(r*r))*dpsidthetaa-(1/r)*dpsidthetadra;

  dvrdthetab = -(1/(r*r))*dpsidtheta2b;
  dvrdrb = (1/(r*r))*dpsidthetab-(1/r)*dpsidthetadrb;

//////////////////////////////////
//  double x = data.r;
//  return sqrt(2.-x)*cos(x);
  if (option == 0) {
    return dvrdthetaa;
  }
  else if (option == 1) {
    return dvrdra;
  }
  else if (option == 2) {
    return dvrdthetab;
  }
  else if (option == 3) {
    return dvrdrb;
  }
  else {
    return 0;
  }
}


static gdouble modified_slip_length(FttCellFace * f, GfsBc * b)
{
  gdouble h = ftt_cell_size (f->cell);
  gdouble Ls0 = gfs_function_face_value (GFS_BC_JOSEPH (b)->lambda, f);
//  gdouble tauc = gfs_function_face_value (GFS_BC_JOSEPH (b)->tauc,f);
  gdouble tauc = gfs_function_face_value (GFS_BC_VALUE (b)->val, f)/Ls0;
//  printf("tauc = %f\n",tauc);

  gdouble theta, beta, U, tau0, tau1, tau2, r01, r02, rc1, rc2, tauavgmod, lambda, dudt, dudn;
  FttVector pos;

  ftt_cell_pos(f->neighbor,&pos);
  
  U = fabs(GFS_VALUE(f->neighbor,b->v)-gfs_function_face_value (GFS_BC_VALUE (b)->val, f));
  dudt = boundary_tangential_gradient(f, b->v);
  dudn = (GFS_VALUE(f->neighbor, b->v)-GFS_VALUE(f->cell,b->v))/h;

  theta = get_dynamic_contact_angle(f->neighbor,b->vofv,0);
  beta = M_PI-theta;
//    tau0 = fabs(dudn); // Simulation stress for Navier BC 

  tau0 = fabs(dudn);//+fabs(dudt); // Simulation stress for Joseph's BC 

// compute new r0 using generalized stokes corner flow solution
  stokes_params data;
  data.phi = theta;
  data.mua = 5.;
  data.mub = 1.;
  data.Lsa = 0.01;
  data.Lsb = 0.01;
  data.alpha = 0.;
  data.U = -U;
  data.r = 1.;
  data.theta = M_PI;

  gdouble dvdna = stokes_v_gradient(data,0); 
  gdouble dvdta = stokes_v_gradient(data,1);

  data.theta = 0.;
  gdouble dvdnb = stokes_v_gradient(data,2);
  gdouble dvdtb = stokes_v_gradient(data,3);
  printf("phi = %f, U = %f, dvdna = %f, dvdta = %f, dvdnb = %f, dvdtb = %f\n",data.phi,data.U,dvdna,dvdta,dvdnb,dvdtb);
//////////////////////////////////////////////
  r01 = (2*sin(theta)*sin(theta))/(theta-sin(theta)*cos(theta))*U/tau0;
  r02 = (2*sin(beta)*sin(beta))/(beta-sin(beta)*cos(beta))*U/tau0;
  rc1 = U/tauc*(2*sin(theta)*sin(theta))/(theta-sin(theta)*cos(theta));
  rc2 = U/tauc*(2*sin(beta)*sin(beta))/(beta-sin(beta)*cos(beta));

  tau1 = tauc*rc1*(1+log(r01/rc1));
  tau2 = tauc*rc2*(1+log(r02/rc2));
  tauavgmod = (tau1+tau2)/(r01+r02);
//  printf("tauavgmod = %f, tau0 = %f\n",tauavgmod,tau0);
  if ((tauavgmod)/tauc > 1 || (tauavgmod) < 0) { 
    lambda = gfs_function_face_value (GFS_BC_VALUE (b)->val, f)/tau0;
//    printf("tauavgmod > tauc: Ls/Ls0 = %f\n",lambda/Ls0);
//    lambda = Ls0/pow(1-0.9999,0.5);
/*    printf("ERROR: (%f,%f) tau_avg > tauc %f > %f, modified LS/Ls0 = %f\n",
            pos.x,pos.y,tauavgmod,tauc,lambda/Ls0); */
  }
  else {
    lambda = Ls0*tauavgmod/tau0;
//    printf("tauavgmod < tauc: Ls/Ls0 = %f\n",lambda/Ls0);
//    lambda = Ls0/pow(1-(tauavgmod)/tauc,0.5);
/*    printf("(%f,%f) Modified Ls/Ls0 = %f\n",pos.x,pos.y,lambda/Ls0); */
  }

  return lambda;
}


static void joseph (FttCellFace * f, GfsBc * b)
{ 

  gdouble h = ftt_cell_size (f->cell);
  gdouble lambda = gfs_function_face_value (GFS_BC_JOSEPH (b)->lambda, f);
  gdouble tauc = gfs_function_face_value (GFS_BC_JOSEPH (b)->tauc,f);


  FttComponent c = f->d % 2;
  FttComponent oc = FTT_ORTHOGONAL_COMPONENT(c);
  FttDirection nd = 2*oc;
  FttCell * nlg = ftt_cell_neighbor (f->cell, nd), * nrg = ftt_cell_neighbor (f->cell, nd + 1);
  FttCell * nl = ftt_cell_neighbor (f->neighbor,nd), * nr = ftt_cell_neighbor (f->neighbor, nd + 1);
  
  guint i;

  FttVector pos, tcp;
  gdouble rc1, rc2, dudn, dudt, tau0, tau1, tau2;
  gdouble Ls0 = lambda, theta, beta;
  gdouble r01, r02, U, tauavgmod, v1, temp_var;

  dudt = boundary_tangential_gradient(f, b->v);
  dudn = (GFS_VALUE(f->neighbor, b->v)-GFS_VALUE(f->cell,b->v))/h;
  v1 = GFS_VALUE(f->neighbor, b->v);

  ftt_cell_pos(f->neighbor,&pos);
  // If cell contains interface, determine modified slip length
  if(b->vofv && is_interfacial(f->neighbor,b->vofv) ) {
/*    gdouble vtop = GFS_STATE (f->neighbor)->f[2].un;
    gdouble vbot = GFS_STATE (f->neighbor)->f[3].un;
    gdouble vcent = GFS_VALUE (f->neighbor,b->v);
    printf("vtop = %f, vbot = %f, vcent = %f\n",vtop,vbot,vcent);
*/
    tcp = find_tcp(f,b);
    if (tcp.x > pos.x + h/2. || tcp.x < pos.x - h/2. || tcp.y > pos.y + h/2. || tcp.y < pos.y - h/2.) {
      /* printf("TCP outside of cell(%f,%f)\n",pos.x,pos.y); */
    }
    else { // TCP inside of cell. Compute modified slip length
      lambda = modified_slip_length(f,b);
    }
  }

  if (!nlg) { // at corner, no left ghost cell
    FttCell * nrrg = ftt_cell_neighbor (nrg, nd+1);
    FttCell * nrr = ftt_cell_neighbor (nr, nd+1);
    
    gdouble vb = (GFS_VALUE(f->cell,b->v)+GFS_VALUE(f->neighbor,b->v))/2.;
    gdouble vr = (GFS_VALUE(nrg,b->v)+GFS_VALUE(nr,b->v))/2.;
    gdouble vrr = (GFS_VALUE(nrrg,b->v)+GFS_VALUE(nrr,b->v))/2.;

    // set ghost cell value 
    GFS_VALUE(f->cell,b->v) =
      (2.*h*(2.*gfs_function_face_value (GFS_BC_VALUE (b)->val, f)-GFS_VALUE(f->neighbor,b->v))+
      lambda*(GFS_VALUE(f->neighbor,b->v)+8.*vr-2.*vrr))/(2.*h+7.*lambda);
  }
  else if (!nrg) { // at corner, no right ghost cell
    FttCell * nllg = ftt_cell_neighbor (nlg, nd);
    FttCell * nll = ftt_cell_neighbor (nl, nd);
    
    gdouble vb = (GFS_VALUE(f->cell,b->v)+GFS_VALUE(f->neighbor,b->v))/2.;
    gdouble vl = (GFS_VALUE(nlg,b->v)+GFS_VALUE(nl,b->v))/2.;
    gdouble vll = (GFS_VALUE(nllg,b->v)+GFS_VALUE(nll,b->v))/2.;
    
    // set ghost cell value 
    GFS_VALUE(f->cell,b->v) = 
      (2.*h*(2.*gfs_function_face_value (GFS_BC_VALUE (b)->val, f)-GFS_VALUE(f->neighbor,b->v))+
      lambda*(7.*GFS_VALUE(f->neighbor,b->v)-8.*vl+2.*vll))/(2.*h+lambda);

  }
  else { // left and right neighbors exist
    gdouble vl = (GFS_VALUE(nlg,b->v)+GFS_VALUE(nl,b->v))/2.;
    gdouble vr = (GFS_VALUE(nrg,b->v)+GFS_VALUE(nr,b->v))/2.;

    // set ghost cell value 
    GFS_VALUE(f->cell,b->v) = 
      (2.*gfs_function_face_value (GFS_BC_VALUE (b)->val, f)*h+
      GFS_VALUE(f->neighbor,b->v)*(2.*lambda-h)+lambda*(vr-vl))/
      (h+2.*lambda);    
  }

}

static void face_joseph (FttCellFace * f, GfsBc * b)
{
  gdouble h = ftt_cell_size (f->cell);
  gdouble lambda = gfs_function_face_value (GFS_BC_JOSEPH (b)->lambda, f);
  FttComponent c = f->d % 2;
  FttComponent oc = FTT_ORTHOGONAL_COMPONENT(c);
  FttDirection nd = 2*oc;
  FttCell * nlg = ftt_cell_neighbor (f->cell, nd), * nrg = ftt_cell_neighbor (f->cell, nd + 1);
  FttCell * nl = ftt_cell_neighbor (f->neighbor,nd), * nr = ftt_cell_neighbor (f->neighbor, nd + 1);
  
  gdouble v1 = GFS_VALUE(f->neighbor,b->v);
  gdouble v0 = gfs_function_face_value (GFS_BC_VALUE (b)->val, f);

  if (!nlg) { // at corner, no left ghost cell
    FttCell * nrrg = ftt_cell_neighbor (nrg, nd+1);
    FttCell * nrr = ftt_cell_neighbor (nr, nd+1);
    
    gdouble vb = (GFS_VALUE(f->cell,b->v)+GFS_VALUE(f->neighbor,b->v))/2.;
    gdouble vr = (GFS_VALUE(nrg,b->v)+GFS_VALUE(nr,b->v))/2.;
    gdouble vrr = (GFS_VALUE(nrrg,b->v)+GFS_VALUE(nrr,b->v))/2.;

    GFS_STATE (f->cell)->f[f->d].v = GFS_STATE (f->neighbor)->f[FTT_OPPOSITE_DIRECTION (f->d)].v = 
      (2.*h*(2.*v0-v1)+lambda*(9.*v1+16.*vr-4.*vrr))/(2.*(2.*h+7.*lambda));
  }
  else if (!nrg) { // at corner, no right ghost cell
    FttCell * nllg = ftt_cell_neighbor (nlg, nd);
    FttCell * nll = ftt_cell_neighbor (nl, nd);
    
    gdouble vb = (GFS_VALUE(f->cell,b->v)+GFS_VALUE(f->neighbor,b->v))/2.;
    gdouble vl = (GFS_VALUE(nlg,b->v)+GFS_VALUE(nl,b->v))/2.;
    gdouble vll = (GFS_VALUE(nllg,b->v)+GFS_VALUE(nll,b->v))/2.;
    
    GFS_STATE (f->cell)->f[f->d].v = GFS_STATE (f->neighbor)->f[FTT_OPPOSITE_DIRECTION (f->d)].v = 
      (2.*h*(4.*v0-v1)+lambda*(15.*v1-16.*vl+4.*vll))/(2.*(2.*h+lambda));
  }
  else { // left and right neighbors exist
    gdouble vl = (GFS_VALUE(nlg,b->v)+GFS_VALUE(nl,b->v))/2.;
    gdouble vr = (GFS_VALUE(nrg,b->v)+GFS_VALUE(nr,b->v))/2.;

  GFS_STATE (f->cell)->f[f->d].v = GFS_STATE (f->neighbor)->f[FTT_OPPOSITE_DIRECTION (f->d)].v = 
    (4.*lambda*v1+2.*v0*h+lambda*(vr-vl))/(2.*h+4.*lambda);
  }
}

static void bc_joseph_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_bc_joseph_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_bc_joseph_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  GfsBcJoseph * bc = GFS_BC_JOSEPH (*o);
  if (bc->lambda == NULL)
    bc->lambda = gfs_function_new (gfs_function_class (), 0.);
  gfs_function_set_units (bc->lambda, 1.);
  gfs_function_read (bc->lambda, gfs_box_domain (GFS_BC (bc)->b->box), fp);
  if (bc->tauc == NULL)
    bc->tauc = gfs_function_new (gfs_function_class (), 10000.);
  gfs_function_set_units (bc->tauc, 1.);
  gfs_function_read (bc->tauc, gfs_box_domain (GFS_BC (bc)->b->box), fp);
}

static void bc_joseph_write (GtsObject * o, FILE * fp)
{  
  (* GTS_OBJECT_CLASS (gfs_bc_joseph_class ())->parent_class->write) (o, fp);
  if (GFS_BC_JOSEPH (o)->lambda)
    gfs_function_write (GFS_BC_JOSEPH (o)->lambda, fp);
  if (GFS_BC_JOSEPH (o)->tauc)
    gfs_function_write (GFS_BC_JOSEPH (o)->tauc, fp);
}

static void gfs_bc_joseph_init (GfsBc * object)
{
  object->bc =                     (FttFaceTraverseFunc) joseph;
  object->homogeneous_bc =         (FttFaceTraverseFunc) homogeneous_dirichlet;
  object->face_bc =                (FttFaceTraverseFunc) face_joseph;
}

static void gfs_bc_joseph_class_init (GtsObjectClass * klass)
{
  klass->read = bc_joseph_read;
  klass->write = bc_joseph_write;
}

GfsBcClass * gfs_bc_joseph_class (void)
{
  static GfsBcClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_bc_joseph_info = {
      "GfsBcJoseph",
      sizeof (GfsBcJoseph),
      sizeof (GfsBcClass),
      (GtsObjectClassInitFunc) gfs_bc_joseph_class_init,
      (GtsObjectInitFunc) gfs_bc_joseph_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_bc_value_class ()),
				  &gfs_bc_joseph_info);
  }

  return klass;
}

/** \endobject{GfsBcJoseph} */

/*----------------------------------------------------------------------------------------------------------*/
/**
 * Navier boundary condition
 * \beginobject{GfsBcNavier}
 */

static void navier (FttCellFace * f, GfsBc * b)
{
  gdouble h = ftt_cell_size (f->cell);
  gdouble lambda = gfs_function_face_value (GFS_BC_NAVIER (b)->lambda, f);

  GFS_VALUE (f->cell, b->v) = 
    (2.*gfs_function_face_value (GFS_BC_VALUE (b)->val, f)*h
     - (h - 2.*lambda)*GFS_VALUE (f->neighbor, b->v))/(h + 2.*lambda);
}

static void face_navier (FttCellFace * f, GfsBc * b)
{
//  printf("face_navier called \n");
  gdouble h = ftt_cell_size (f->cell);
  gdouble lambda = gfs_function_face_value (GFS_BC_NAVIER (b)->lambda, f);
  GFS_STATE (f->cell)->f[f->d].v = GFS_STATE (f->neighbor)->f[FTT_OPPOSITE_DIRECTION (f->d)].v = 
    (gfs_function_face_value (GFS_BC_VALUE (b)->val, f)*h + 
     2.*lambda*GFS_VALUE (f->neighbor, b->v))/(h + 2.*lambda);
}

static void bc_navier_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_bc_navier_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_bc_navier_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  GfsBcNavier * bc = GFS_BC_NAVIER (*o);
  if (bc->lambda == NULL)
    bc->lambda = gfs_function_new (gfs_function_class (), 0.);
  gfs_function_set_units (bc->lambda, 1.);
  gfs_function_read (bc->lambda, gfs_box_domain (GFS_BC (bc)->b->box), fp);
}

static void bc_navier_write (GtsObject * o, FILE * fp)
{  
  (* GTS_OBJECT_CLASS (gfs_bc_navier_class ())->parent_class->write) (o, fp);
  if (GFS_BC_NAVIER (o)->lambda)
    gfs_function_write (GFS_BC_NAVIER (o)->lambda, fp);
}

static void gfs_bc_navier_init (GfsBc * object)
{
  object->bc =                     (FttFaceTraverseFunc) navier;
  object->homogeneous_bc =         (FttFaceTraverseFunc) homogeneous_dirichlet;
  object->face_bc =                (FttFaceTraverseFunc) face_navier;
}

static void gfs_bc_navier_class_init (GtsObjectClass * klass)
{
  klass->read = bc_navier_read;
  klass->write = bc_navier_write;
}

GfsBcClass * gfs_bc_navier_class (void)
{
  static GfsBcClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_bc_navier_info = {
      "GfsBcNavier",
      sizeof (GfsBcNavier),
      sizeof (GfsBcClass),
      (GtsObjectClassInitFunc) gfs_bc_navier_class_init,
      (GtsObjectInitFunc) gfs_bc_navier_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_bc_value_class ()),
          &gfs_bc_navier_info);
  }

  return klass;
}

/** \endpbject{GfsBcNavier} */
/*----------------------------------------------------------------------------------------------------------*/



/**
 * One of the boundaries of a #GfsBox.
 * \beginobject{GfsBoundary}
 */

static void insert_bc (GfsVariable * v, GtsObject * o, GHashTable * unique)
{
  g_hash_table_insert (unique, o, o);
}

static void gfs_boundary_destroy (GtsObject * object)
{
  GfsBoundary * boundary = GFS_BOUNDARY (object);
  GfsDomain * domain = gfs_box_domain (boundary->box);

  if (domain)
    gfs_domain_forget_boundary (domain, boundary);
  if (boundary->root) {
    if (domain == NULL) /* domain has been destroyed */
      ftt_cell_destroy (boundary->root, NULL, NULL);
    else
      ftt_cell_destroy (boundary->root, (FttCellCleanupFunc) gfs_cell_cleanup, domain);
  }
  boundary->box->neighbor[FTT_OPPOSITE_DIRECTION (boundary->d)] = NULL;

  gts_object_destroy (GTS_OBJECT (boundary->default_bc));
  if (boundary->bc) {
    GHashTable * unique = g_hash_table_new (NULL, NULL);
    /* make sure that the pointers are unique before destroying them */
    g_hash_table_foreach (boundary->bc, (GHFunc) insert_bc, unique);
    g_hash_table_foreach (unique, (GHFunc) gts_object_destroy, NULL);
    g_hash_table_destroy (unique);
    g_hash_table_destroy (boundary->bc);
  }

  (* GTS_OBJECT_CLASS (gfs_boundary_class ())->parent_class->destroy) (object);
}

static void match (FttCell * cell, GfsBoundary * boundary)
{
  FttCell * neighbor = ftt_cell_neighbor (cell, boundary->d);
  FttCell * parent = ftt_cell_parent (cell);
  guint level = ftt_cell_level (cell);

  cell->flags |= GFS_FLAG_BOUNDARY;
  if (parent && GFS_CELL_IS_GRADIENT_BOUNDARY (parent))
    cell->flags |= GFS_FLAG_GRADIENT_BOUNDARY;
  if (neighbor == NULL || ftt_cell_level (neighbor) < level) {
    if (FTT_CELL_IS_ROOT (cell)) {
      g_assert (cell == boundary->root);
      boundary->root = NULL;
    }
    ftt_cell_destroy (cell, (FttCellCleanupFunc) gfs_cell_cleanup, gfs_box_domain (boundary->box));
    boundary->changed = TRUE;
    return;
  }
  if (ftt_cell_level (neighbor) == level) {
    GfsSolidVector * s = GFS_STATE (neighbor)->solid;

    if (s && s->s[FTT_OPPOSITE_DIRECTION (boundary->d)] == 0.) {
      if (FTT_CELL_IS_ROOT (cell)) {
	g_assert (cell == boundary->root);
	boundary->root = NULL;
      }
      ftt_cell_destroy (cell, (FttCellCleanupFunc) gfs_cell_cleanup,
			gfs_box_domain (boundary->box));
      boundary->changed = TRUE;
      return;
    }
    if (s) {
      FttDirection d;
      FttComponent c;
      GfsSolidVector * t;

      if (GFS_STATE (cell)->solid == NULL)
	GFS_STATE (cell)->solid = g_malloc0 (sizeof (GfsSolidVector));
      t = GFS_STATE (cell)->solid;
      t->a = s->a;
      for (d = 0; d < FTT_NEIGHBORS; d++)
	if (d/2 == boundary->d/2)
	  t->s[d] = s->s[FTT_OPPOSITE_DIRECTION (d)];
	else
	  t->s[d] = s->s[d];
      for (c = 0; c < FTT_DIMENSION; c++)
	if (c == boundary->d/2) {
	  FttVector p1, p2;
	  ftt_cell_pos (cell, &p1);
	  ftt_cell_pos (neighbor, &p2);
	  (&t->cm.x)[c] = (&p1.x)[c] + (&p2.x)[c] - (&s->cm.x)[c];
	  (&t->ca.x)[c] = (&p1.x)[c] + (&p2.x)[c] - (&s->ca.x)[c];
	}
	else {
	  (&t->cm.x)[c] = (&s->cm.x)[c];
	  (&t->ca.x)[c] = (&s->ca.x)[c];
	}
    }
    else if (GFS_STATE (cell)->solid != NULL) {
      g_free (GFS_STATE (cell)->solid);
      GFS_STATE (cell)->solid = NULL;
    }      
    if (FTT_CELL_IS_LEAF (cell) && !FTT_CELL_IS_LEAF (neighbor)) {
      GfsDomain * domain = gfs_box_domain (boundary->box);
      ftt_cell_refine_single (cell, domain->cell_init, domain->cell_init_data);
      boundary->changed = TRUE;
    }
  }
  else
    g_assert_not_reached ();
  if (!FTT_CELL_IS_LEAF (cell))
    level++;
  if (level > boundary->depth)
    boundary->depth = level;
}

static void boundary_match (GfsBoundary * boundary)
{
  if (boundary->root == NULL) {
    GfsBox * box = boundary->box;
    GfsDomain * domain = gfs_box_domain (box);
    boundary->root = ftt_cell_new ((FttCellInitFunc) gfs_cell_init, domain);
    FTT_ROOT_CELL (boundary->root)->parent = box;
    ftt_cell_set_level (boundary->root, ftt_cell_level (box->root));
    ftt_cell_set_neighbor_match (boundary->root, box->root, boundary->d, 
				 (FttCellInitFunc) gfs_cell_init, domain);
    FttVector pos;
    ftt_cell_pos (box->root, &pos);
    gdouble size = ftt_cell_size (box->root);
    FttDirection d = FTT_OPPOSITE_DIRECTION (boundary->d);
    pos.x += rpos[d].x*size;
    pos.y += rpos[d].y*size;
    pos.z += rpos[d].z*size;
    ftt_cell_set_pos (boundary->root, &pos);
  }

  guint l = ftt_cell_level (boundary->root);
  
  boundary->changed = FALSE;
  boundary->depth = l;
  while (boundary->root && l <= boundary->depth) {
    ftt_cell_traverse_boundary (boundary->root, boundary->d,
				FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL, l,
				(FttCellTraverseFunc) match, boundary);
    l++;
  }
  if (boundary->root && boundary->changed)
    ftt_cell_flatten (boundary->root, boundary->d, (FttCellCleanupFunc) gfs_cell_cleanup, 
		      gfs_box_domain (boundary->box));
}

static void is_extra (GfsVariable * v, GfsBc * bc, gboolean * extra)
{
  if (bc->extra)
    *extra = TRUE;
}

static void write_extra (GfsVariable * v, GfsBc * bc, FILE * fp)
{
  if (bc->extra) {
    if(GTS_OBJECT (bc)->klass->write) {
      (* GTS_OBJECT (bc)->klass->write) (GTS_OBJECT (bc), fp);
      fputc ('\n', fp);
    }
  }
}

static void gfs_boundary_write (GtsObject * o, FILE * fp)
{
  GfsBoundary * b = GFS_BOUNDARY (o);
  gboolean any_extra = FALSE;

  if (GTS_OBJECT_CLASS (gfs_boundary_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_boundary_class ())->parent_class->write) 
      (o, fp);

  g_hash_table_foreach (b->bc, (GHFunc) is_extra, &any_extra);
  if (any_extra) {
    fputs (" {\n", fp);
    g_hash_table_foreach (b->bc, (GHFunc) write_extra, fp);
    fputc ('}', fp);
  }
}

static gboolean boundary_read_extra_bc (GfsBoundary * b, GtsFile * fp)
{
  gboolean ret = FALSE;

  if (fp->type != '{')
    return ret;

  fp->scope_max++;
  gts_file_next_token (fp);

  while (fp->type != GTS_ERROR && fp->type != '}') {
    if (fp->type == '\n') {
      gts_file_next_token (fp);
      continue;
    }
    if (fp->type != GTS_STRING) {
      gts_file_error (fp, "expecting a keyword");
      return ret;
    }
    else {
      GtsObjectClass * klass = gfs_object_class_from_name (fp->token->str);
      GtsObject * object;
      
      if (klass == NULL) {
	gts_file_error (fp, "unknown keyword `%s'", fp->token->str);
	return ret;
      }
      else if (!gts_object_class_is_from_class (klass, gfs_bc_class ())) {
	gts_file_error (fp, "`%s' is not a GfsBc", fp->token->str);
	return ret;
      }

      object = gts_object_new (klass);
      g_assert (klass->read);
      GFS_BC (object)->b = b;
      GFS_BC (object)->extra = TRUE;
      (* klass->read) (&object, fp);
      if (fp->type == GTS_ERROR) {
	gts_object_destroy (object);
	return ret;
      }

      gfs_boundary_add_bc (b, GFS_BC (object));
      ret = TRUE;
    }
  }
  if (fp->type != '}') {
    gts_file_error (fp, "expecting a closing brace");
    return ret;
  }
  fp->scope_max--;
  gts_file_next_token (fp);

  return ret;
}

static void gfs_boundary_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_boundary_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_boundary_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsBoundary * b = GFS_BOUNDARY (*o);
  GfsVariable ** v = gfs_domain_velocity (gfs_box_domain (b->box));
  if (v)
    gfs_boundary_add_bc (b, gfs_bc_value_new (gfs_bc_dirichlet_class (), v[b->d/2], NULL, FALSE));

  boundary_read_extra_bc (GFS_BOUNDARY (*o), fp);
}

static void gfs_boundary_class_init (GfsBoundaryClass * klass)
{
  klass->match = boundary_match;

  GTS_OBJECT_CLASS (klass)->write =   gfs_boundary_write;
  GTS_OBJECT_CLASS (klass)->read =    gfs_boundary_read;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_boundary_destroy;
}

static void gfs_boundary_init (GfsBoundary * b)
{
  b->type = GFS_BOUNDARY_CENTER_VARIABLE;
  b->bc = g_hash_table_new (g_str_hash, g_str_equal);
  gfs_boundary_set_default_bc (b, gfs_bc_new (gfs_bc_class (), NULL, FALSE));
}

GfsBoundaryClass * gfs_boundary_class (void)
{
  static GfsBoundaryClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_boundary_info = {
      "GfsBoundary",
      sizeof (GfsBoundary),
      sizeof (GfsBoundaryClass),
      (GtsObjectClassInitFunc) gfs_boundary_class_init,
      (GtsObjectInitFunc) gfs_boundary_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (gts_object_class (), &gfs_boundary_info);
  }

  return klass;
}

/**
 * gfs_boundary_new:
 * @klass: a #GfsBoundaryClass.
 * @box: a #GfsBox.
 * @d: a direction.
 *
 * Creates a new boundary of type @klass for @box in direction @d.
 *
 * This function fails if @box has already a boundary in direction @d.
 *
 * Returns: a new #GfsBoundary.
 */
GfsBoundary * gfs_boundary_new (GfsBoundaryClass * klass,
				GfsBox * box,
				FttDirection d)
{
  GfsBoundary * boundary;

  g_return_val_if_fail (box != NULL, NULL);
  g_return_val_if_fail (d < FTT_NEIGHBORS, NULL);
  g_return_val_if_fail (box->neighbor[d] == NULL, NULL);

  boundary = GFS_BOUNDARY (gts_object_new (GTS_OBJECT_CLASS (klass)));
  boundary->box = box;
  box->neighbor[d] = GTS_OBJECT (boundary);
  boundary->d = FTT_OPPOSITE_DIRECTION (d);
  if (box->root)
    boundary_match (boundary);

  return boundary;
}

/**
 * gfs_boundary_update:
 * @boundary: a #GfsBoundary.
 *
 * Calls the @update() method of @boundary.
 */
void gfs_boundary_update (GfsBoundary * boundary)
{
  g_return_if_fail (boundary != NULL);

  if (GFS_BOUNDARY_CLASS (GTS_OBJECT (boundary)->klass)->update)
    (* GFS_BOUNDARY_CLASS (GTS_OBJECT (boundary)->klass)->update) (boundary);
}

/**
 * gfs_boundary_send:
 * @boundary: a #GfsBoundary.
 *
 * Calls the @send() method of @boundary.
 */
void gfs_boundary_send (GfsBoundary * boundary)
{
  g_return_if_fail (boundary != NULL);

  if (GFS_BOUNDARY_CLASS (GTS_OBJECT (boundary)->klass)->send)
    (* GFS_BOUNDARY_CLASS (GTS_OBJECT (boundary)->klass)->send) (boundary);
}

/**
 * gfs_boundary_receive:
 * @boundary: a #GfsBoundary.
 * @flags: the traversal flags.
 * @max_depth: the maximum depth of the traversal.
 *
 * Calls the @receive() method of @boundary.
 */
void gfs_boundary_receive (GfsBoundary * boundary,
			   FttTraverseFlags flags,
			   gint max_depth)
{
  g_return_if_fail (boundary != NULL);

  if (GFS_BOUNDARY_CLASS (GTS_OBJECT (boundary)->klass)->receive)
    (* GFS_BOUNDARY_CLASS (GTS_OBJECT (boundary)->klass)->receive)
      (boundary, flags, max_depth);
}

/**
 * gfs_boundary_synchronize:
 * @boundary: a #GfsBoundary.
 *
 * Calls the @synchronize() method of @boundary.
 */
void gfs_boundary_synchronize (GfsBoundary * boundary)
{
  g_return_if_fail (boundary != NULL);

  if (GFS_BOUNDARY_CLASS (GTS_OBJECT (boundary)->klass)->synchronize)
    (* GFS_BOUNDARY_CLASS (GTS_OBJECT (boundary)->klass)->synchronize) (boundary);
}

/**
 * gfs_boundary_lookup:
 * @b: a #GfsBoundary.
 * @v: a #GfsVariable.
 *
 * Returns: the #GfsBc associated with @b and @v.
 */
GfsBc * gfs_boundary_lookup_bc (GfsBoundary * b, GfsVariable * v)
{
  GfsBc * bv;

  g_return_val_if_fail (b != NULL, NULL);
  g_return_val_if_fail (v != NULL, NULL);

  if (!v->name || !(bv = g_hash_table_lookup (b->bc, v->name))) {
    if (v->default_bc) {
      bv = v->default_bc;
      bv->b = b;
    }
    else
      bv = b->default_bc;
    bv->v = v;
  }
  return bv;
}

/**
 * gfs_boundary_set_default_bc:
 * @b: a #GfsBoundary.
 * @bc: a #GfsBc.
 *
 * Sets the default boundary condition for @b to @bc.
 */
void gfs_boundary_set_default_bc (GfsBoundary * b, GfsBc * bc)
{
  g_return_if_fail (b != NULL);
  g_return_if_fail (bc != NULL);
  g_return_if_fail (bc->b == NULL || bc->b == b);

  if (b->default_bc)
    gts_object_destroy (GTS_OBJECT (b->default_bc));
  b->default_bc = bc;
  bc->b = b;
}

/**
 * gfs_variable_set_default_bc:
 * @v: a #GfVariable.
 * @bc: a #GfsBc.
 *
 * Sets the default boundary condition for @v to @bc.
 */
void gfs_variable_set_default_bc (GfsVariable * v, GfsBc * bc)
{
  g_return_if_fail (v != NULL);
  g_return_if_fail (bc != NULL);
  g_return_if_fail (bc->v == NULL || bc->v == v);

  if (v->default_bc)
    gts_object_destroy (GTS_OBJECT (v->default_bc));
  v->default_bc = bc;
  bc->v = v;
}

/**
 * gfs_boundary_add_bc:
 * @b: a #GfsBoundary.
 * @bc: a #GfsBc.
 *
 * Adds boundary condition @bc to @b.
 */
void gfs_boundary_add_bc (GfsBoundary * b, GfsBc * bc)
{
  GfsBc * old;
  
  g_return_if_fail (b != NULL);
  g_return_if_fail (bc != NULL);
  g_return_if_fail (bc->v != NULL);
  g_return_if_fail (bc->v->name != NULL);
  g_return_if_fail (bc->b == NULL || bc->b == b);

  old = g_hash_table_lookup (b->bc, bc->v->name);
  if (!old || !old->extra) {
    if (old) gts_object_destroy (GTS_OBJECT (old));
    g_hash_table_insert (b->bc, bc->v->name, bc);
    bc->b = b;
  }
  else
    gts_object_destroy (GTS_OBJECT (bc));
}

/** \endobject{GfsBoundary} */

/**
 * Constant inflow boundary.
 * \beginobject{GfsBoundaryInflowConstant}
 */

static GtsColor inflow_color (GtsObject * o)
{
  GtsColor c = { 0., 0., 1. }; /* blue */

  return c;
}

static void inflow_constant_write (GtsObject * o, FILE * fp)
{
  if (GTS_OBJECT_CLASS (gfs_boundary_inflow_constant_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_boundary_inflow_constant_class ())->parent_class->write) 
      (o, fp);

  gfs_function_write (GFS_BOUNDARY_INFLOW_CONSTANT (o)->un, fp);
}

static void inflow_constant_read (GtsObject ** o, GtsFile * fp)
{
  GfsBoundary * b = GFS_BOUNDARY (*o);  
  FttComponent c;
  GfsFunction * un = GFS_BOUNDARY_INFLOW_CONSTANT (*o)->un;
  GfsVariable ** v;

  if (GTS_OBJECT_CLASS (gfs_boundary_inflow_constant_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_boundary_inflow_constant_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  gfs_function_read (un, gfs_box_domain (b->box), fp);
  gfs_function_set_units (un, 1.);

  v = gfs_domain_velocity (gfs_box_domain (b->box));
  for (c = 0; c < FTT_DIMENSION; c++)
    if (c == b->d/2)
      gfs_boundary_add_bc (b, gfs_bc_value_new (gfs_bc_dirichlet_class (),
						v[c], un, FALSE));
    else
      gfs_boundary_add_bc (b, gfs_bc_value_new (gfs_bc_dirichlet_class (),
						v[c], NULL, FALSE));
}

static void gfs_boundary_inflow_constant_class_init (GtsObjectClass * klass)
{
  klass->read    = inflow_constant_read;
  klass->write   = inflow_constant_write;
  klass->color   = inflow_color;
}

static void gfs_boundary_inflow_constant_init (GfsBoundaryInflowConstant * object)
{
  object->un = gfs_function_new (gfs_function_class (), 0.);
}

GfsBoundaryInflowConstantClass * gfs_boundary_inflow_constant_class (void)
{
  static GfsBoundaryInflowConstantClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_boundary_inflow_constant_info = {
      "GfsBoundaryInflowConstant",
      sizeof (GfsBoundaryInflowConstant),
      sizeof (GfsBoundaryInflowConstantClass),
      (GtsObjectClassInitFunc) gfs_boundary_inflow_constant_class_init,
      (GtsObjectInitFunc) gfs_boundary_inflow_constant_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_boundary_class ()),
				  &gfs_boundary_inflow_constant_info);
  }

  return klass;
}

/** \endobject{GfsBoundaryInflowConstant} */

/**
 * Outflow boundary.
 * \beginobject{GfsBoundaryOutflow}
 */

static GtsColor outflow_color (GtsObject * o)
{
  GtsColor c = { 0., 1., 0. }; /* green */

  return c;
}

static void outflow_read (GtsObject ** o, GtsFile * fp)
{
  GfsBoundary * b = GFS_BOUNDARY (*o);
  GfsDomain * domain;
  GfsVariable ** v;

  if (GTS_OBJECT_CLASS (gfs_boundary_outflow_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_boundary_outflow_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  domain = gfs_box_domain (b->box);
  v = gfs_domain_velocity (domain);
  gfs_boundary_add_bc (b, gfs_bc_value_new (gfs_bc_neumann_class (),
					    v[b->d/2],
					    NULL, FALSE));
  gfs_boundary_add_bc (b, gfs_bc_value_new (gfs_bc_dirichlet_class (),
					    gfs_variable_from_name (domain->variables, "P"),
					    NULL, FALSE));
}

static void gfs_boundary_outflow_class_init (GfsBoundaryClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read   = outflow_read;
  GTS_OBJECT_CLASS (klass)->color  = outflow_color;
}

GfsBoundaryOutflowClass * gfs_boundary_outflow_class (void)
{
  static GfsBoundaryOutflowClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_boundary_outflow_info = {
      "GfsBoundaryOutflow",
      sizeof (GfsBoundary),
      sizeof (GfsBoundaryOutflowClass),
      (GtsObjectClassInitFunc) gfs_boundary_outflow_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_boundary_class ()),
				  &gfs_boundary_outflow_info);
  }

  return klass;
}

/** \endobject{GfsBoundaryOutflow} */

/**
 * 
 * \beginobject{GfsBoundaryGradient}
 */

static GtsColor gradient_color (GtsObject * o)
{
  GtsColor c = { 1., 1., 0. }; /* red-green */

  return c;
}

static void set_gradient_boundary (FttCell * cell)
{
  cell->flags |= GFS_FLAG_GRADIENT_BOUNDARY;
}

static void gradient_read (GtsObject ** o, GtsFile * fp)
{
  GfsBoundary * b = GFS_BOUNDARY (*o);

  if (GTS_OBJECT_CLASS (gfs_boundary_gradient_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_boundary_gradient_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  ftt_cell_traverse (b->root, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
		     (FttCellTraverseFunc) set_gradient_boundary, NULL);
}

static void gfs_boundary_gradient_class_init (GfsBoundaryClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read   = gradient_read;
  GTS_OBJECT_CLASS (klass)->color  = gradient_color;
}

GfsBoundaryClass * gfs_boundary_gradient_class (void)
{
  static GfsBoundaryClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_boundary_gradient_info = {
      "GfsBoundaryGradient",
      sizeof (GfsBoundary),
      sizeof (GfsBoundaryClass),
      (GtsObjectClassInitFunc) gfs_boundary_gradient_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_boundary_class ()),
				  &gfs_boundary_gradient_info);
  }

  return klass;
}

/** \endobject{GfsBoundaryGradient} */

/**
 * Periodic boundary.
 * \beginobject{GfsBoundaryPeriodic}
 */

static void boundary_periodic_destroy (GtsObject * object)
{
  GfsBoundaryPeriodic * boundary = GFS_BOUNDARY_PERIODIC (object);

  g_array_free (boundary->sndbuf, TRUE);
  g_array_free (boundary->rcvbuf, TRUE);
  
  (* GTS_OBJECT_CLASS (gfs_boundary_periodic_class ())->parent_class->destroy) 
    (object);
}

static void boundary_periodic_read (GtsObject ** object, GtsFile * fp)
{
  boundary_periodic_destroy (*object);
}

static void center_periodic (FttCellFace * face, GfsBc * b)
{
  GfsBoundaryPeriodic * boundary_periodic = GFS_BOUNDARY_PERIODIC (b->b);

  g_assert (boundary_periodic->sndcount < boundary_periodic->sndbuf->len);
  g_assert (ftt_face_type (face) == FTT_FINE_FINE);
  g_assert (!FTT_CELL_IS_LEAF (face->cell) || FTT_CELL_IS_LEAF (face->neighbor));
  g_array_index (boundary_periodic->sndbuf, gdouble, boundary_periodic->sndcount++) =
    GFS_VALUE (face->neighbor, b->v);
}

static void face_periodic (FttCellFace * face, GfsBc * b)
{
  GfsBoundaryPeriodic * boundary_periodic = GFS_BOUNDARY_PERIODIC (b->b);

  g_assert (boundary_periodic->sndcount < boundary_periodic->sndbuf->len);
  g_array_index (boundary_periodic->sndbuf, gdouble, boundary_periodic->sndcount++) =
    GFS_STATE (face->neighbor)->f[FTT_OPPOSITE_DIRECTION (face->d)].v;
}

static void boundary_size (FttCell * cell, guint * count)
{
  (*count)++;
}

static void set_buffers_size (GfsBoundaryPeriodic * boundary)
{
  guint count = 0;

  ftt_cell_traverse (GFS_BOUNDARY (boundary)->root, 
		     FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
		     (FttCellTraverseFunc) boundary_size, &count);
  g_array_set_size (boundary->rcvbuf, count);
  g_array_set_size (boundary->sndbuf, count);
}

static void boundary_tree (FttCell * cell, GfsBoundaryPeriodic * boundary)
{
  gdouble is_leaf = FTT_CELL_IS_LEAF (cell);

  if (boundary->sndcount == boundary->sndbuf->len)
    g_array_append_val (boundary->sndbuf, is_leaf);
  else
    g_array_index (boundary->sndbuf, gdouble, boundary->sndcount) = is_leaf;
  boundary->sndcount++;

  if (!is_leaf) {
    FttCellChildren child;
    guint i, n;

    n = ftt_cell_children_direction (cell, GFS_BOUNDARY (boundary)->d, &child);
    for (i = 0; i < n; i++) {
      /* fixme: using a gdouble to store (and MPI transfer) a boolean is wasteful... */
      gdouble is_destroyed = (child.c[i] == NULL);
      
      if (boundary->sndcount == boundary->sndbuf->len)
	g_array_append_val (boundary->sndbuf, is_destroyed);
      else
	g_array_index (boundary->sndbuf, gdouble, boundary->sndcount) = is_destroyed;
      boundary->sndcount++;
    }
    for (i = 0; i < n; i++)
      if (child.c[i])
	boundary_tree (child.c[i], boundary);
  }
}

static void periodic_match (GfsBoundary * boundary)
{
  (* gfs_boundary_class ()->match) (boundary);

  g_assert (GFS_BOUNDARY_PERIODIC (boundary)->sndcount == 0);
  if (boundary->root)
    boundary_tree (boundary->root, GFS_BOUNDARY_PERIODIC (boundary));
}

static void send (GfsBoundary * bb)
{
  GfsBoundaryPeriodic * boundary = GFS_BOUNDARY_PERIODIC (bb);
  g_assert (boundary->matching);
  GfsBoundaryPeriodic * matching = 
    GFS_BOUNDARY_PERIODIC (boundary->matching->neighbor[boundary->d]);

  g_assert (GFS_IS_BOUNDARY_PERIODIC (matching));
  g_assert (boundary->sndcount <= boundary->sndbuf->len);
  
  if (GFS_BOUNDARY (boundary)->type == GFS_BOUNDARY_MATCH_VARIABLE) {
    if (boundary->sndcount > matching->rcvbuf->len)
      g_array_set_size (matching->rcvbuf, boundary->sndcount);
  }
  memcpy (matching->rcvbuf->data, boundary->sndbuf->data, boundary->sndcount*sizeof (gdouble));
}

static void center_update (FttCell * cell,
			   GfsBoundaryPeriodic * boundary)
{
  g_assert (boundary->rcvcount < boundary->rcvbuf->len);
  GFS_VALUE (cell, GFS_BOUNDARY (boundary)->v) =
    g_array_index (boundary->rcvbuf, gdouble, boundary->rcvcount++);
}

static void face_update (FttCellFace * face,
			 GfsBoundaryPeriodic * boundary)
{
  g_assert (boundary->rcvcount < boundary->rcvbuf->len);
  GFS_STATE (face->cell)->f[face->d].v = 
    g_array_index (boundary->rcvbuf, gdouble, boundary->rcvcount++);
}

static void match_ignore (GfsBoundaryPeriodic * boundary)
{
  gboolean is_leaf;

  g_assert (boundary->rcvcount < boundary->rcvbuf->len);
  is_leaf = g_array_index (boundary->rcvbuf, gdouble, boundary->rcvcount++);

  if (!is_leaf) {
    gboolean is_destroyed[FTT_CELLS/2];
    guint i;

    for (i = 0; i < FTT_CELLS/2; i++) {
      g_assert (boundary->rcvcount < boundary->rcvbuf->len);
      is_destroyed[i] = g_array_index (boundary->rcvbuf, gdouble, boundary->rcvcount++);
    }
    for (i = 0; i < FTT_CELLS/2; i++)
      if (!is_destroyed[i])
	match_ignore (boundary);
  }
}

static void match_update (FttCell * cell,
			  GfsBoundaryPeriodic * boundary)
{
  gboolean is_leaf;

  g_assert (boundary->rcvcount < boundary->rcvbuf->len);
  is_leaf = g_array_index (boundary->rcvbuf, gdouble, boundary->rcvcount++);

  if (!is_leaf) {
    GfsDomain * domain = gfs_box_domain (GFS_BOUNDARY (boundary)->box);
    FttCellChildren child;
    gboolean is_destroyed[FTT_CELLS/2];
    guint i, n;

    if (FTT_CELL_IS_LEAF (cell)) {
      FttCell * neighbor = ftt_cell_neighbor (cell, GFS_BOUNDARY (boundary)->d);

      g_assert (neighbor);
      ftt_cell_refine_single (cell, domain->cell_init, domain->cell_init_data);
      if (FTT_CELL_IS_LEAF (neighbor))
	ftt_cell_refine_single (neighbor, domain->cell_init, domain->cell_init_data);
      /* what about solid fractions? */
      GFS_BOUNDARY (boundary)->changed = TRUE;
    }
    n = ftt_cell_children_direction (cell, GFS_BOUNDARY (boundary)->d, &child);
    for (i = 0; i < n; i++) {
      g_assert (boundary->rcvcount < boundary->rcvbuf->len);
      is_destroyed[i] = g_array_index (boundary->rcvbuf, gdouble, boundary->rcvcount++);
      if (is_destroyed[i] && child.c[i]) {
	ftt_cell_destroy (child.c[i], (FttCellCleanupFunc) gfs_cell_cleanup, domain);
	child.c[i] = NULL;
	GFS_BOUNDARY (boundary)->changed = TRUE;
      }
    }
    for (i = 0; i < n; i++)
      if (!is_destroyed[i]) {
	if (child.c[i])
	  match_update (child.c[i], boundary);
	else
	  match_ignore (boundary);
      }
  }
}

static void receive (GfsBoundary * bb,
		     FttTraverseFlags flags,
		     gint max_depth)
{
  GfsBoundaryPeriodic * boundary = GFS_BOUNDARY_PERIODIC (bb);

  boundary->rcvcount = 0;
  switch (GFS_BOUNDARY (boundary)->type) {
  case GFS_BOUNDARY_FACE_VARIABLE:
    ftt_face_traverse_boundary (GFS_BOUNDARY (boundary)->root,
				GFS_BOUNDARY (boundary)->d,
				FTT_PRE_ORDER, flags, max_depth,
				(FttFaceTraverseFunc) face_update, boundary);
    break;

  case GFS_BOUNDARY_MATCH_VARIABLE:
    match_update (GFS_BOUNDARY (boundary)->root, boundary);
    ftt_cell_flatten (GFS_BOUNDARY (boundary)->root, 
		      GFS_BOUNDARY (boundary)->d,
		      (FttCellCleanupFunc) gfs_cell_cleanup, 
		      gfs_box_domain (GFS_BOUNDARY (boundary)->box));
    break;

  default:
    ftt_cell_traverse (GFS_BOUNDARY (boundary)->root,
		       FTT_PRE_ORDER, flags, max_depth,
		       (FttCellTraverseFunc) center_update, boundary);
  }
}

static void synchronize (GfsBoundary * bb)
{
  GfsBoundaryPeriodic * boundary = GFS_BOUNDARY_PERIODIC (bb);

  boundary->sndcount = 0;
  if (bb->type == GFS_BOUNDARY_MATCH_VARIABLE)
    set_buffers_size (boundary);
}

static GtsColor periodic_color (GtsObject * o)
{
  GtsColor c = { 1., 0., 0. }; /* red */

  return c;
}

static void gfs_boundary_periodic_class_init (GfsBoundaryClass * klass)
{
  GfsBoundaryClass * parent_class = GFS_BOUNDARY_CLASS (klass);

  parent_class->match             = periodic_match;
  parent_class->send              = send;
  parent_class->receive           = receive;
  parent_class->synchronize       = synchronize;

  GTS_OBJECT_CLASS (klass)->color =   periodic_color;
  GTS_OBJECT_CLASS (klass)->destroy = boundary_periodic_destroy;
  GTS_OBJECT_CLASS (klass)->read    = boundary_periodic_read;
}

static void gfs_boundary_periodic_init (GfsBoundaryPeriodic * boundary)
{
  GfsBc * b = GFS_BOUNDARY (boundary)->default_bc;

  b->bc                = (FttFaceTraverseFunc) center_periodic;
  b->homogeneous_bc    = (FttFaceTraverseFunc) center_periodic;
  b->face_bc           = (FttFaceTraverseFunc) face_periodic;

  boundary->sndbuf = g_array_new (FALSE, FALSE, sizeof (gdouble));
  boundary->rcvbuf = g_array_new (FALSE, FALSE, sizeof (gdouble));
  boundary->sndcount = boundary->rcvcount = 0;

  boundary->rotate = 0.;
}

GfsBoundaryClass * gfs_boundary_periodic_class (void)
{
  static GfsBoundaryClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_boundary_periodic_info = {
      "GfsBoundaryPeriodic",
      sizeof (GfsBoundaryPeriodic),
      sizeof (GfsBoundaryClass),
      (GtsObjectClassInitFunc) gfs_boundary_periodic_class_init,
      (GtsObjectInitFunc) gfs_boundary_periodic_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_boundary_class ()),
				  &gfs_boundary_periodic_info);
  }

  return klass;
}

/**
 * gfs_boundary_periodic_new:
 * @klass: a #GfsBoundaryClass.
 * @box: a #GfsBox.
 * @d: a #FttDirection.
 * @matching: a #GfsBox.
 *
 * Returns: a new #GfsBoundaryPeriodic connecting @box with @matching in direction @d.
 */
GfsBoundaryPeriodic * gfs_boundary_periodic_new (GfsBoundaryClass * klass,
						 GfsBox * box,
						 FttDirection d,
						 GfsBox * matching)
{
  GfsBoundaryPeriodic * boundary;

  boundary = GFS_BOUNDARY_PERIODIC (gfs_boundary_new (klass, box, d));
  set_buffers_size (boundary);
  boundary->matching = matching;
  boundary->d = FTT_OPPOSITE_DIRECTION (d);

  return boundary;
}

static void center_periodic_rotate (FttCellFace * face, GfsBc * b)
{
  GfsBoundaryPeriodic * boundary_periodic = GFS_BOUNDARY_PERIODIC (b->b);

  g_assert (boundary_periodic->sndcount < boundary_periodic->sndbuf->len);
  g_assert (ftt_face_type (face) == FTT_FINE_FINE);
  g_assert (!FTT_CELL_IS_LEAF (face->cell) || FTT_CELL_IS_LEAF (face->neighbor));

  if (b->v->component < 2) { /* 2D-vector-rotation only */
    FttComponent c = FTT_ORTHOGONAL_COMPONENT (b->v->component);
    g_assert (b->v->vector[c]);
    g_array_index (boundary_periodic->sndbuf, gdouble, boundary_periodic->sndcount++) =
      b->v->orientation*(b->v->even ? 1. : boundary_periodic->rotate)*
      GFS_VALUE (face->neighbor, b->v->vector[c]);
  }
  else
    g_array_index (boundary_periodic->sndbuf, gdouble, boundary_periodic->sndcount++) =
      GFS_VALUE (face->neighbor, b->v);
}

static void face_periodic_rotate (FttCellFace * face, GfsBc * b)
{
  GfsBoundaryPeriodic * boundary_periodic = GFS_BOUNDARY_PERIODIC (b->b);

  g_assert (boundary_periodic->sndcount < boundary_periodic->sndbuf->len);
  FttDirection d = FTT_OPPOSITE_DIRECTION (face->d);
  if (b->v->component < 2) { /* 2D-vector-rotation only */
    FttComponent c = FTT_ORTHOGONAL_COMPONENT (b->v->component);
    g_assert (d < 4);
    g_assert (b->v->domain->has_rotated_bc);
    g_assert (b->v->face[c][d]);
    g_array_index (boundary_periodic->sndbuf, gdouble, boundary_periodic->sndcount++) =
      (2.*c - 1.)*boundary_periodic->rotate*GFS_VALUE (face->neighbor, b->v->face[c][d]);
  }
  else
    g_array_index (boundary_periodic->sndbuf, gdouble, boundary_periodic->sndcount++) =
      GFS_STATE (face->neighbor)->f[d].v;
}

/**
 * gfs_boundary_periodic_rotate:
 * @boundary: a #GfsBoundaryPeriodic.
 * @rotate: a #FttDirection.
 * @orientation: the orientation (+1 or -1).
 *
 * Rotates #boundary according to @rotate and @orientation.
 */
void gfs_boundary_periodic_rotate (GfsBoundaryPeriodic * boundary,
				   FttDirection rotate,
				   gdouble orientation)
{
  g_return_if_fail (boundary != NULL);

  boundary->d = rotate;
  boundary->rotate = orientation;
  gfs_box_domain (GFS_BOUNDARY (boundary)->box)->has_rotated_bc = TRUE;

  GfsBc * b = GFS_BOUNDARY (boundary)->default_bc;
  b->bc = b->homogeneous_bc = (FttFaceTraverseFunc) center_periodic_rotate;
  b->face_bc = (FttFaceTraverseFunc) face_periodic_rotate;
}

/**
 * gfs_boundary_periodic_rotate_new:
 * @klass: a #GfsBoundaryClass.
 * @box: a #GfsBox.
 * @d: a #FttDirection.
 * @matching: a #GfsBox.
 * @rotate: a #FttDirection.
 * @orientation: the orientation (+1 or -1).
 *
 * Returns: a new "rotated" #GfsBoundaryPeriodic connecting @box in
 * direction @d with @matching in direction @rotate, oriented using
 * @orientation.
 */
GfsBoundaryPeriodic * gfs_boundary_periodic_rotate_new (GfsBoundaryClass * klass,
							GfsBox * box,
							FttDirection d,
							GfsBox * matching,
							FttDirection rotate,
							gdouble orientation)
{
  GfsBoundaryPeriodic * boundary;

  boundary = gfs_boundary_periodic_new (klass, box, d, matching);
  gfs_boundary_periodic_rotate (boundary, rotate, orientation);

  return boundary;
}

/** \endobject{GfsBoundaryPeriodic} */

/**
 * Link between #GfsBox.
 * \beginobject{GfsGEdge}
 */

static void gfs_gedge_write (GtsObject * object, FILE * fp)
{
  fprintf (fp, " %s", ftt_direction_name [GFS_GEDGE (object)->d]);
  if (GFS_GEDGE (object)->rotate < FTT_NEIGHBORS)
    fprintf (fp, " %s", ftt_direction_name [GFS_GEDGE (object)->rotate]);
}

static void gfs_gedge_read (GtsObject ** o, GtsFile * fp)
{
  GfsGEdge * e = GFS_GEDGE (*o);

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (direction)");
    return;
  }
  e->d = ftt_direction_from_name (fp->token->str);
  if (e->d >= FTT_NEIGHBORS) {
    gts_file_error (fp, "unknown direction `%s'", fp->token->str);
    e->d = 0;
    return;
  }
  gts_file_next_token (fp);
  if (fp->type == GTS_STRING) {
    e->rotate = ftt_direction_from_name (fp->token->str);
    if (e->rotate >= FTT_NEIGHBORS) {
      gts_file_error (fp, "unknown direction `%s'", fp->token->str);
      e->rotate = -1;
      return;
    }
    gts_file_next_token (fp);
  }
}

static void gfs_gedge_class_init (GtsObjectClass * klass)
{
  klass->write = gfs_gedge_write;
  klass->read = gfs_gedge_read;
}

static void gfs_gedge_init (GfsGEdge * object)
{
  object->d = object->rotate = FTT_NEIGHBORS;
}

GfsGEdgeClass * gfs_gedge_class (void)
{
  static GfsGEdgeClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_gedge_info = {
      "GfsGEdge",
      sizeof (GfsGEdge),
      sizeof (GfsGEdgeClass),
      (GtsObjectClassInitFunc) gfs_gedge_class_init,
      (GtsObjectInitFunc) gfs_gedge_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_gedge_class ()),
				  &gfs_gedge_info);
  }

  return klass;
}

/**
 * gfs_gedge_link_boxes:
 * @edge: a #GfsGEdge.
 *
 * Links the two boxes connected by @edge. The boxes are set as their
 * respective neighbors in the direction defined by @edge (note that
 * their relative positions are not set).
 */
void gfs_gedge_link_boxes (GfsGEdge * edge)
{
  GfsBox * b1, * b2;

  g_return_if_fail (edge != NULL);
  g_return_if_fail (GTS_GEDGE (edge)->n1 != NULL);
  g_return_if_fail (GTS_GEDGE (edge)->n2 != NULL);
  g_return_if_fail (edge->d >= 0 && edge->d < FTT_NEIGHBORS);

  b1 = GFS_BOX (GTS_GEDGE (edge)->n1);
  b2 = GFS_BOX (GTS_GEDGE (edge)->n2);

  g_return_if_fail (b1->neighbor[edge->d] == NULL);
  
  if (edge->rotate < FTT_NEIGHBORS) {
    g_return_if_fail (b2->neighbor[edge->rotate] == NULL);
    gfs_boundary_periodic_rotate_new (gfs_boundary_periodic_class (),
				      b1, edge->d, b2, edge->rotate,  1.);
    gfs_boundary_periodic_rotate_new (gfs_boundary_periodic_class (),
				      b2, edge->rotate, b1, edge->d, -1.);
  }
  else {
    g_return_if_fail (b2->neighbor[FTT_OPPOSITE_DIRECTION (edge->d)] == NULL);
    
    GtsObject * periodic = GTS_OBJECT (b1);
    gdouble * p1 = &FTT_ROOT_CELL (b1->root)->pos.x;
    gdouble * p2 = &FTT_ROOT_CELL (b2->root)->pos.x;
    gdouble sign = edge->d % 2 ? 1. : -1.;
    FttComponent c = edge->d/2;
    if (p1[c] != G_MAXDOUBLE && p2[c] != G_MAXDOUBLE && sign*(p2[c] - p1[c]) > 0.)
      periodic = GTS_OBJECT (b2);
    else
      while (periodic && GFS_IS_BOX (periodic) && GFS_BOX (periodic) != b2)
	periodic = GFS_BOX (periodic)->neighbor[FTT_OPPOSITE_DIRECTION (edge->d)];
    
    if (GFS_BOX (periodic) == b2) {
      gfs_boundary_periodic_new (gfs_boundary_periodic_class (), b1, edge->d, b2);
      gfs_boundary_periodic_new (gfs_boundary_periodic_class (), b2, 
				 FTT_OPPOSITE_DIRECTION (edge->d), b1);
    }
    else {
      GfsDomain * domain = gfs_box_domain (b1);
      if (domain->pid < 0 || b1->pid == b2->pid)
	ftt_cell_set_neighbor (b1->root, b2->root, edge->d, 
			       (FttCellInitFunc) gfs_cell_init, domain);
      b1->neighbor[edge->d] = GTS_OBJECT (b2);
      b2->neighbor[FTT_OPPOSITE_DIRECTION (edge->d)] = GTS_OBJECT (b1);
    }
  }
}

/**
 * gfs_gedge_new:
 * @klass: a #GfsGEdgeClass.
 * @b1: a #GfsBox.
 * @b2: another #GfsBox.
 * @d: a direction.
 *
 * Returns: a new #GfsGEdge linking @b1 to @b2 in direction @d. The
 * boxes then need to be linked using gfs_gedge_link_boxes().
 */
GfsGEdge * gfs_gedge_new (GfsGEdgeClass * klass,
			  GfsBox * b1, GfsBox * b2,
			  FttDirection d)
{
  GfsGEdge * edge;

  g_return_val_if_fail (klass != NULL, NULL);
  g_return_val_if_fail (b1 != NULL, NULL);
  g_return_val_if_fail (b2 != NULL, NULL);
  g_return_val_if_fail (d >= 0 && d < FTT_NEIGHBORS, NULL);

  edge = GFS_GEDGE (gts_gedge_new (GTS_GEDGE_CLASS (klass),
				   GTS_GNODE (b1), GTS_GNODE (b2)));
  edge->d = d;

  return edge;
}

/** \endobject{GfsGEdge} */

/**
 * The box making up a #GfsDomain.
 * \beginobject{GfsBox}
 */

static void gfs_box_destroy (GtsObject * object)
{
  GfsBox * box = GFS_BOX (object);
  FttDirection d;

  if (box->root) {
    GfsDomain * domain = gfs_box_domain (box);
    if (domain == NULL) /* domain has been destroyed */
      ftt_cell_destroy (box->root, NULL, NULL);
    else
      ftt_cell_destroy (box->root, (FttCellCleanupFunc) gfs_cell_cleanup, domain);
  }
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY (box->neighbor[d]))
      gts_object_destroy (box->neighbor[d]);
    else if (GFS_IS_BOX (box->neighbor[d])) {
      g_assert (GFS_BOX (box->neighbor[d])->neighbor[FTT_OPPOSITE_DIRECTION (d)] == 
		GTS_OBJECT (box));
      GFS_BOX (box->neighbor[d])->neighbor[FTT_OPPOSITE_DIRECTION (d)] = NULL;
    }

  (* GTS_OBJECT_CLASS (gfs_box_class ())->parent_class->destroy) (object);
}

static void box_size (FttCell * cell, guint * size)
{
  (*size)++;
}

static void gfs_box_write (GtsObject * object, FILE * fp)
{
  GfsBox * box = GFS_BOX (object);
  FttVector pos;
  FttDirection d;
  guint size = 0;
  GfsDomain * domain = gfs_box_domain (box);

  ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		     (FttCellTraverseFunc) box_size, &size);
  ftt_cell_pos (box->root, &pos);
  fprintf (fp, "%s { id = %u pid = %d size = %u x = %g y = %g z = %g",
	   object->klass->info.name, box->id, box->pid, size, pos.x, pos.y, pos.z);
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY (box->neighbor[d])) {
      fprintf (fp, " %s = %s",
	       ftt_direction_name[d],
	       box->neighbor[d]->klass->info.name);
      if (box->neighbor[d]->klass->write)
	(* box->neighbor[d]->klass->write) (box->neighbor[d], fp);
    }
  fputs (" }", fp);
  if (domain != NULL && domain->max_depth_write > -2) {
    fputs (" {\n", fp);
    if (domain->binary)
      ftt_cell_write_binary (box->root, domain->max_depth_write, fp, 
			     (FttCellWriteFunc) gfs_cell_write_binary, domain->variables_io);
    else
      ftt_cell_write (box->root, domain->max_depth_write, fp, 
		      (FttCellWriteFunc) gfs_cell_write, domain->variables_io);
    fputc ('}', fp);
  }
}

static void gfs_box_read (GtsObject ** o, GtsFile * fp)
{
  GfsBox * b = GFS_BOX (*o);
  GtsObjectClass * klass;
  gboolean class_changed = FALSE;
  FttVector pos = {0., 0., 0.};
  GtsFileVariable var[] = {
    {GTS_UINT,   "id",     TRUE, &b->id},
    {GTS_INT,    "pid",    TRUE, &b->pid},
    {GTS_UINT,   "size",   TRUE, &b->size},
    {GTS_DOUBLE, "x",      TRUE, &pos.x},
    {GTS_DOUBLE, "y",      TRUE, &pos.y},
    {GTS_DOUBLE, "z",      TRUE, &pos.z},
    {GTS_FILE,   "right",  TRUE},
    {GTS_FILE,   "left",   TRUE},
    {GTS_FILE,   "top",    TRUE},
    {GTS_FILE,   "bottom", TRUE},
#if (!FTT_2D)
    {GTS_FILE,   "front",  TRUE},
    {GTS_FILE,   "back",   TRUE},
#endif /* 3D */
    {GTS_NONE}
  };
  GtsFileVariable * v;
  gfloat weight;
  GfsDomain * domain = GTS_OBJECT (*o)->reserved;

  if (domain == NULL) {
    g_assert (GTS_SLIST_CONTAINEE (b)->containers &&
	      !GTS_SLIST_CONTAINEE (b)->containers->next);
    domain = GFS_DOMAIN (GTS_SLIST_CONTAINEE (b)->containers->data);
  }
  else
    gts_container_add (GTS_CONTAINER (domain), GTS_CONTAINEE (b));

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (GfsBoxClass)");
    return;
  }
  klass = gfs_object_class_from_name (fp->token->str);
  if (klass == NULL) {
    gts_file_error (fp, "unknown class `%s'", fp->token->str);
    return;
  }
  if (!gts_object_class_is_from_class (klass, gfs_box_class ())) {
    gts_file_error (fp, "`%s' is not a GfsBox", fp->token->str);
    return;
  }
  if (klass != (*o)->klass) {
    *o = gts_object_new (klass);
    gts_object_destroy (GTS_OBJECT (b));
    b = GFS_BOX (*o);
    gts_container_add (GTS_CONTAINER (domain), GTS_CONTAINEE (b));
    class_changed = TRUE;
  }
  gts_file_next_token (fp);

  g_assert (b->root == NULL);
  b->root = ftt_cell_new ((FttCellInitFunc) gfs_cell_init, domain);

  weight = gts_gnode_weight (GTS_GNODE (b));
  gts_file_assign_start (fp, var);
  if (fp->type == GTS_ERROR)
    return;

  while ((v = gts_file_assign_next (fp, var)))
    if (v->type == GTS_FILE) {
      GtsObjectClass * boundary_class = gfs_object_class_from_name (fp->token->str);
      GtsObject * boundary;
	
      if (boundary_class == NULL) {
	gts_file_error (fp, "unknown class `%s'", fp->token->str);
	return;
      }
      if (!gts_object_class_is_from_class (boundary_class, gfs_boundary_class ())) {
	gts_file_error (fp, "`%s' is not a GfsBoundary", fp->token->str);
	return;
      }
      boundary = GTS_OBJECT (gfs_boundary_new (GFS_BOUNDARY_CLASS (boundary_class),
					       b, ftt_direction_from_name (v->name)));
      gts_file_next_token (fp);
      if (boundary_class->read)
	(* boundary_class->read) (&boundary, fp);
    }
  
  if (fp->type == '{') {
    FttCell * root;

    fp->scope_max++;
    if (domain->binary) {
      if (gts_file_getc (fp) != '\n') {
      	gts_file_error (fp, "expecting a newline");
      	return;
      }
      root = ftt_cell_read_binary (fp, (FttCellReadFunc) gfs_cell_read_binary, domain);
      if (fp->type == GTS_ERROR)
	return;
      gts_file_next_token (fp);
    }
    else {
      gts_file_first_token_after (fp, '\n');
      root = ftt_cell_read (fp, (FttCellReadFunc) gfs_cell_read, domain);
    }
    fp->scope_max--;

    if (domain->pid >= 0 && b->pid != domain->pid)
      /* ignore data of boxes belonging to other PEs */
      ftt_cell_destroy (root, (FttCellCleanupFunc) gfs_cell_cleanup, domain);
    else {
      ftt_cell_destroy (b->root, (FttCellCleanupFunc) gfs_cell_cleanup, domain);
      b->root = root;
      FttDirection d;
      for (d = 0; d < FTT_NEIGHBORS; d++)
	if (GFS_IS_BOUNDARY (b->neighbor[d])) {
	  GfsBoundary * boundary = GFS_BOUNDARY (b->neighbor[d]);
	  
	  ftt_cell_set_neighbor_match (boundary->root, b->root, boundary->d, 
				       (FttCellInitFunc) gfs_cell_init, domain);
	}
    }
	
    if (fp->type == GTS_ERROR)
      return;
    if (fp->type != '}') {
      gts_file_error (fp, "expecting a closing brace");
      return;
    }
    gts_file_next_token (fp);
  }

  FTT_ROOT_CELL (b->root)->parent = b;
  if (ftt_cell_level (b->root) != domain->rootlevel) {
    FttDirection d;

    ftt_cell_set_level (b->root, domain->rootlevel);
    for (d = 0; d < FTT_NEIGHBORS; d++)
      if (GFS_IS_BOUNDARY (b->neighbor[d]))
	ftt_cell_set_level (GFS_BOUNDARY (b->neighbor[d])->root, domain->rootlevel);
  }

  if (var[3].set || var[4].set || var[5].set) {
    gdouble size = ftt_cell_size (b->root);
    FttDirection d;
    ftt_cell_set_pos (b->root, &pos);
    for (d = 0; d < FTT_NEIGHBORS; d++)
      if (GFS_IS_BOUNDARY (b->neighbor[d])) {
	FttVector bpos = pos;
	bpos.x += rpos[d].x*size;
	bpos.y += rpos[d].y*size;
	bpos.z += rpos[d].z*size;
	ftt_cell_set_pos (GFS_BOUNDARY (b->neighbor[d])->root, &bpos);
      }
  }
  else /* position is not set */
    FTT_ROOT_CELL (b->root)->pos.x = G_MAXDOUBLE;

  /* updates weight of domain */
  GTS_WGRAPH (domain)->weight += gts_gnode_weight (GTS_GNODE (b)) - weight;

  if (class_changed && klass->read)
    (* klass->read) (o, fp);
}

static gfloat gfs_box_weight (GtsGNode * node)
{
  GfsBox * box = GFS_BOX (node);

  if (box->size >= 0)
    return box->size;
  else {
    guint size = 0;

    if (box->root)
      ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			 (FttCellTraverseFunc) box_size, &size);
    return size;
  }
}

static void gfs_box_class_init (GfsBoxClass * klass)
{
  GTS_GNODE_CLASS (klass)->weight = gfs_box_weight;

  GTS_OBJECT_CLASS (klass)->destroy = gfs_box_destroy;
  GTS_OBJECT_CLASS (klass)->write = gfs_box_write;
  GTS_OBJECT_CLASS (klass)->read = gfs_box_read;
}

static void gfs_box_init (GfsBox * box)
{
  static guint id = 1;

  box->id = id++;
  box->pid = -1;
  box->size = -1;
}

GfsBoxClass * gfs_box_class (void)
{
  static GfsBoxClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_box_info = {
      "GfsBox",
      sizeof (GfsBox),
      sizeof (GfsBoxClass),
      (GtsObjectClassInitFunc) gfs_box_class_init,
      (GtsObjectInitFunc) gfs_box_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_gnode_class ()),
				  &gfs_box_info);
  }

  return klass;
}

GfsBox * gfs_box_new (GfsBoxClass * klass)
{
  GfsBox * object;

  object = GFS_BOX (gts_object_new (GTS_OBJECT_CLASS (klass)));

  return object;
}

/** \endobject{GfsBox} */
