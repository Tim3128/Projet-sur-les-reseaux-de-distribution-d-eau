///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  MONITEUR D'ENCHAINEMENT POUR LE CALCUL DE L'EQUILIBRE D'UN RESEAU D'EAU  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// --------------------------------------
// Dimensionnement de l'espace de travail
// --------------------------------------

   stacksize(10000000);

// ------------------------------------------
// Fonctions fournies dans le cadre du projet
// ------------------------------------------

   // Donnees du problemes

   exec('Probleme_R.sce');
   exec('Structures_R.sce');
   
   // Affichage des resultats

   exec('Visualg.sci');
   
   // Verification  des resultats

   exec('HydrauliqueP.sci');
   exec('HydrauliqueD.sci');
   exec('Verification.sci');

// ------------------------------------------
// Fonctions a ecrire dans le cadre du projet
// ------------------------------------------

   // ---> Charger les fonctions  associees a l'oracle du probleme,
   //      aux algorithmes d'optimisation et de recherche lineaire.
   //
   // Exemple : la fonction "optim" de Scilab
   //
   exec('OraclePG.sci');
   exec('OraclePH.sci');
   exec('OracleDG.sci');
   exec('OracleDH.sci');
   exec('Gradient_F.sci');
   exec('Gradient_V.sci');
   exec('Gradient_F.sci');
   exec('Polack_Ribiere.sci');
   exec('BFGS.sci');
   exec('Newton.sci');
   exec('Optim_Scilab.sci');
   titrgr = "Fonction Newton à pas variable sur le probleme primal";

// ------------------------------
// Initialisation de l'algorithme
// ------------------------------

   // La dimension (n-md) est celle du probleme primal

   xini = 0.1 * rand(n-md,1);
   lambdaini = 0.1 * rand(m-mr,1);

// ----------------------------
// Minimisation proprement dite
// ----------------------------

   // Exemple : la fonction "optim" de Scilab
   //[fopt, xopt, gopt] = optim(OraclePG, xini);
   //[fopt,xopt,gopt] = Gradient_F(OraclePH,xini);
   //[fopt,xopt,gopt] = Gradient_V(OraclePH,xini);
   //[fopt,xopt,gopt] = Polack_Ribiere(OraclePH,xini);
   //[fopt,xopt,gopt] = BFGS(OraclePH,xini);
   [fopt,xopt,gopt] = Newton(OraclePH,xini);

   // Minimisation du duale
   
   //[fopt,lambdaopt,gopt] = Gradient_V(OracleDG,lambdaini);
   //[fopt,lamdaopt,gopt] = Polack_Ribiere(OracleDH,lambdaini);
   //[fopt,lambdaopt,gopt] = BFGS(OracleDH,lambdaini);
   // [fopt,lambdaopt,gopt] = Newton(OracleDH,lambdaini);

// --------------------------
// Verification des resultats
// --------------------------

   [q,z,f,p] = HydrauliqueP(xopt);
   
   Verification(q,z,f,p);
   
   //[q,z,f,p] = HydrauliqueD(lambdaopt);
   
   //Verification(q,z,f,p);



//
