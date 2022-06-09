// Copyright 2019-2022 Cambridge Quantum Computing
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <catch2/catch_test_macros.hpp>
#include <cstdlib>

#include "../testutil.hpp"
#include "Utils/Constants.hpp"
#include "Utils/CosSinDecomposition.hpp"
#include "Utils/MatrixAnalysis.hpp"

namespace tket {
namespace test_CosSinDecomposition {

static bool cs_matrices_ok(const Eigen::MatrixXd &c, const Eigen::MatrixXd &s) {
  unsigned n = c.rows();
  if (!c.isDiagonal()) return false;
  if (!s.isDiagonal()) return false;
  if (!Eigen::MatrixXcd::Identity(n, n).isApprox(c * c + s * s)) return false;
  for (unsigned i = 0; i < n; i++) {
    if (c(i, i) < 0) return false;
    if (s(i, i) < 0) return false;
  }
  for (unsigned i = 0; i + 1 < n; i++) {
    if (c(i, i) > c(i + 1, i + 1) + ERR_EPS) return false;
  }
  return true;
}

static void test_csd(const Eigen::MatrixXcd &U) {
  auto [l0, l1, r0, r1, c, s] = CS_decomp(U);

  unsigned N = U.rows();
  unsigned n = N / 2;

  Eigen::MatrixXcd L = Eigen::MatrixXcd::Zero(N, N);
  L.topLeftCorner(n, n) = l0;
  L.bottomRightCorner(n, n) = l1;
  Eigen::MatrixXd CS = Eigen::MatrixXd::Zero(N, N);
  CS.topLeftCorner(n, n) = c;
  CS.topRightCorner(n, n) = -s;
  CS.bottomLeftCorner(n, n) = s;
  CS.bottomRightCorner(n, n) = c;
  Eigen::MatrixXcd R = Eigen::MatrixXcd::Zero(N, N);
  R.topLeftCorner(n, n) = r0;
  R.bottomRightCorner(n, n) = r1;

  CHECK(U.isApprox(L * CS * R));
  CHECK(cs_matrices_ok(c, s));
  CHECK(is_unitary(l0));
  CHECK(is_unitary(l1));
  CHECK(is_unitary(r0));
  CHECK(is_unitary(r1));
}

SCENARIO("CosSinDecomposition") {
  GIVEN("Fixed 2x2 unitary") {
    Eigen::Matrix2cd U2 = Eigen::Matrix2cd::Zero();
    U2(0, 0) = {0.2817184155378645, 0.3796799045050548};
    U2(0, 1) = {0.7710111478974819, -0.4266376850205267};
    U2(1, 0) = {0.7782225612421542, 0.41333708959585147};
    U2(1, 1) = {-0.2751580580853349, 0.38446084145051745};
    test_csd(U2);
  }
  GIVEN("Fixed 4x4 unitary") {
    Eigen::Matrix4cd U4 = Eigen::Matrix4cd::Zero();
    U4(0, 0) = {0.13679196211550004, 0.7041818405777518};
    U4(0, 1) = {0.16128284952412175, 0.4788107685040558};
    U4(0, 2) = {0.16249995838978795, 0.4402567889298545};
    U4(0, 3) = {-0.006489742649654795, 0.09934581762655259};
    U4(1, 0) = {0.34888657967038395, -0.280031360791438};
    U4(1, 1) = {-0.6171950395501815, -0.11399476635854978};
    U4(1, 2) = {0.21117623780276323, 0.5743024177112399};
    U4(1, 3) = {-0.07590232512624884, 0.16048749027768602};
    U4(2, 0) = {-0.11707259941153635, 0.2149549505736465};
    U4(2, 1) = {-0.17223822670989994, -0.20714710122894883};
    U4(2, 2) = {0.4097974958734698, 0.014760960771141208};
    U4(2, 3) = {-0.11817260277711523, -0.8278863507565566};
    U4(3, 0) = {-0.07045672403903053, 0.4694685417883851};
    U4(3, 1) = {-0.36803970710776407, -0.3778550152052355};
    U4(3, 2) = {0.3031834158741162, -0.38115219489704427};
    U4(3, 3) = {0.31296551789335747, 0.40157552846117056};
    test_csd(U4);
  }
  GIVEN("Fixed 8x8 unitary") {
    Eigen::MatrixXcd U8 = Eigen::MatrixXcd::Zero(8, 8);
    U8(0, 0) = {0.13316830729343884, -0.19488348823413731};
    U8(0, 1) = {-0.31512556534985303, -0.3540942467438625};
    U8(0, 2) = {0.32045718213700497, 0.10718780060303738};
    U8(0, 3) = {-0.15560380327827494, 0.2945106546484882};
    U8(0, 4) = {-0.10736181286021396, -0.3770716188213386};
    U8(0, 5) = {0.32708787433308645, 0.04264199120153849};
    U8(0, 6) = {0.42321768751249467, -0.18355814821549069};
    U8(0, 7) = {0.12903764590617614, -0.0499676951979884};
    U8(1, 0) = {-0.1894601584306871, -0.07007929300249177};
    U8(1, 1) = {0.35716640739552763, 0.08571620831845965};
    U8(1, 2) = {-0.05610454132434403, -0.04798510507074549};
    U8(1, 3) = {0.2120948031489579, 0.28002281918740934};
    U8(1, 4) = {-0.5285096169885108, -0.4524925178620848};
    U8(1, 5) = {-0.20509425512538984, 0.15032847199758032};
    U8(1, 6) = {0.09980158567192993, 0.28545060715797066};
    U8(1, 7) = {-0.15000140591829803, 0.18098157671159207};
    U8(2, 0) = {0.25471935188898676, 0.1520758154092516};
    U8(2, 1) = {-0.5596750195621968, -0.22584557737663913};
    U8(2, 2) = {-0.27683501138691635, 0.2896354328407763};
    U8(2, 3) = {0.3260585632028611, 0.14805416892608741};
    U8(2, 4) = {-0.07309997411917406, -0.07551215157679726};
    U8(2, 5) = {-0.07071951860033429, -0.07862231787876223};
    U8(2, 6) = {-0.326131579921589, 0.31511354734055863};
    U8(2, 7) = {-0.0076616743680552315, 0.17618923892243882};
    U8(3, 0) = {-0.4270750048857905, 0.422945967115775};
    U8(3, 1) = {-0.05464735369810664, 0.009111430554056325};
    U8(3, 2) = {0.2787619854618884, 0.23404713940582147};
    U8(3, 3) = {0.025835239080968693, -0.26595745106472896};
    U8(3, 4) = {-0.029012032328860876, -0.1636886170847321};
    U8(3, 5) = {0.42037387602904686, -0.15812760365096107};
    U8(3, 6) = {-0.0571543881170567, 0.1854237509979565};
    U8(3, 7) = {-0.38809593205096926, -0.11893663966454963};
    U8(4, 0) = {-0.17613994315390888, -0.2430614486685144};
    U8(4, 1) = {0.16114823805916656, -0.2690764363901607};
    U8(4, 2) = {0.3515132298101505, -0.11131296696777837};
    U8(4, 3) = {0.1314261306284719, -0.24871643877108404};
    U8(4, 4) = {0.06245010497269367, -0.15564556659819814};
    U8(4, 5) = {0.030797896123874003, -0.2896285932317556};
    U8(4, 6) = {-0.3539805917205374, 0.12719220883389995};
    U8(4, 7) = {0.49362735032749333, 0.31358137467517955};
    U8(5, 0) = {-0.20536807458986261, -0.0424855362407279};
    U8(5, 1) = {0.09882112625229451, -0.1819724858784419};
    U8(5, 2) = {0.10109888342811012, 0.5514704742015339};
    U8(5, 3) = {-0.01234038384108932, 0.10628791738700805};
    U8(5, 4) = {0.001569069423803214, 0.10943644845734934};
    U8(5, 5) = {-0.5546576910750309, -0.4305701476190374};
    U8(5, 6) = {0.15879830567593642, -0.18989483221207706};
    U8(5, 7) = {-0.05108058787980293, -0.1358251735615456};
    U8(6, 0) = {-0.05290044361649104, 0.48445699202511755};
    U8(6, 1) = {0.11695864942757844, 0.34074446455202195};
    U8(6, 2) = {0.05199157333705276, 0.2652187261390794};
    U8(6, 3) = {-0.10595776033487594, 0.3825026935850879};
    U8(6, 4) = {-0.04574000023916003, -0.019048940801146677};
    U8(6, 5) = {0.11489677584634068, 0.08506093727879346};
    U8(6, 6) = {-0.15351256113945266, -0.10641384772090345};
    U8(6, 7) = {0.5863103748411339, -0.02441075987325043};
    U8(7, 0) = {0.04239722400064411, 0.2842802197407317};
    U8(7, 1) = {0.008317021052256901, 0.044503429884961046};
    U8(7, 2) = {-0.23601175044094616, -0.09119553843999911};
    U8(7, 3) = {0.4662501874862328, -0.31703659236298676};
    U8(7, 4) = {0.2656553762562674, -0.45869929569315404};
    U8(7, 5) = {-0.11243770688739829, -0.051804608049269396};
    U8(7, 6) = {0.31146474700382726, -0.3432003414714381};
    U8(7, 7) = {0.14938070011657345, 0.0016608482084349248};
    test_csd(U8);
  }
  GIVEN("Some special matrices") {
    Eigen::Matrix4cd U4;
    // clang-format off
    U4 << 0, 0, 1, 0,
          0, 0, 0, 1,
          1, 0, 0, 0,
          0, 1, 0, 0;
    test_csd(U4);
    U4 << 0, 0, 0, 1,
          0, 0, 1, 0,
          0, 1, 0, 0,
          1, 0, 0, 0;
    test_csd(U4);
    // clang-format on
    Eigen::MatrixXcd U8 = Eigen::MatrixXcd::Zero(8, 8);
    for (unsigned i = 0; i < 4; i++) {
      U8(i, 4 + i) = 1;
      U8(4 + i, i) = 1;
    }
    test_csd(U8);
  }
  GIVEN("Random unitaries") {
    for (unsigned n = 2; n <= 8; n += 2) {
      for (unsigned i = 0; i < 100; i++) {
        Eigen::MatrixXcd U = random_unitary(n, 100 * n + i);
        test_csd(U);
      }
    }
  }
  GIVEN("Direct sums") {
    for (unsigned n = 1; n <= 4; n++) {
      for (unsigned i = 0; i < 10; i++) {
        Eigen::MatrixXcd U = Eigen::MatrixXcd::Zero(2 * n, 2 * n);
        U.topLeftCorner(n, n) = random_unitary(n, 10 * n + i);
        U.bottomRightCorner(n, n) = random_unitary(n, 100 + 10 * n + i);
        test_csd(U);
      }
    }
  }
  GIVEN("Kronecker products") {
    for (unsigned n = 1; n <= 4; n++) {
      for (unsigned i = 0; i < 10; i++) {
        Eigen::MatrixXcd U = random_unitary(n, 10 * n + i);
        Eigen::MatrixXcd V = random_unitary(2, 100 + 10 * n + i);
        test_csd(Eigen::kroneckerProduct(U, V));
        test_csd(Eigen::kroneckerProduct(V, U));
      }
    }
  }
  GIVEN("Identity") {
    for (unsigned n = 2; n <= 8; n += 2) {
      test_csd(Eigen::MatrixXcd::Identity(n, n));
    }
  }
  GIVEN("Direct sum of random with identity") {
    for (unsigned n = 1; n <= 4; n++) {
      const Eigen::MatrixXcd I_2n = Eigen::MatrixXcd::Identity(2 * n, 2 * n);
      for (unsigned i = 0; i < 10; i++) {
        Eigen::MatrixXcd U = I_2n;
        const Eigen::MatrixXcd R = random_unitary(n, 10 * n + i);
        U.topLeftCorner(n, n) = R;
        test_csd(U);
        Eigen::MatrixXcd V = I_2n;
        V.bottomRightCorner(n, n) = R;
        test_csd(V);
      }
    }
  }
  GIVEN("Kronecker product of random with identity") {
    for (unsigned n = 1; n <= 4; n++) {
      const Eigen::MatrixXcd I_n = Eigen::MatrixXcd::Identity(n, n);
      for (unsigned i = 0; i < 10; i++) {
        Eigen::MatrixXcd R = random_unitary(2, 10 * n + i);
        test_csd(Eigen::kroneckerProduct(I_n, R));
        test_csd(Eigen::kroneckerProduct(R, I_n));
      }
    }
  }
  GIVEN("Kronecker product of random with X") {
    for (unsigned n = 1; n <= 4; n++) {
      Eigen::Matrix2cd X = Eigen::Matrix2cd::Zero();
      X(0, 1) = 1;
      X(1, 0) = 1;
      for (unsigned i = 0; i < 10; i++) {
        Eigen::MatrixXcd R = random_unitary(n, 10 * n + i);
        test_csd(Eigen::kroneckerProduct(X, R));
        test_csd(Eigen::kroneckerProduct(R, X));
      }
    }
  }
}

}  // namespace test_CosSinDecomposition
}  // namespace tket
