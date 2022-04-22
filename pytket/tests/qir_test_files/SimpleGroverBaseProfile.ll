; ModuleID = 'qat-link'
source_filename = "qat-link"

%Qubit = type opaque
%Result = type opaque

define void @Microsoft__Quantum__Samples__SimpleGrover__SearchForMarkedInput__Interop() local_unnamed_addr #0 {
entry:
  tail call void @__quantum__qis__h__body(%Qubit* null)
  tail call void @__quantum__qis__h__body(%Qubit* nonnull inttoptr (i64 1 to %Qubit*))
  tail call void @__quantum__qis__x__body(%Qubit* nonnull inttoptr (i64 2 to %Qubit*))
  tail call void @__quantum__qis__h__body(%Qubit* nonnull inttoptr (i64 2 to %Qubit*))
  tail call void @__quantum__qis__x__body(%Qubit* null)
  tail call void @__quantum__qis__h__body(%Qubit* nonnull inttoptr (i64 2 to %Qubit*))
  tail call void @__quantum__qis__t__adj(%Qubit* null)
  tail call void @__quantum__qis__t__adj(%Qubit* nonnull inttoptr (i64 1 to %Qubit*))
  tail call void @__quantum__qis__cnot__body(%Qubit* nonnull inttoptr (i64 2 to %Qubit*), %Qubit* null)
  tail call void @__quantum__qis__t__body(%Qubit* null)
  tail call void @__quantum__qis__cnot__body(%Qubit* nonnull inttoptr (i64 1 to %Qubit*), %Qubit* nonnull inttoptr (i64 2 to %Qubit*))
  tail call void @__quantum__qis__cnot__body(%Qubit* nonnull inttoptr (i64 1 to %Qubit*), %Qubit* null)
  tail call void @__quantum__qis__t__body(%Qubit* nonnull inttoptr (i64 2 to %Qubit*))
  tail call void @__quantum__qis__t__adj(%Qubit* null)
  tail call void @__quantum__qis__cnot__body(%Qubit* nonnull inttoptr (i64 1 to %Qubit*), %Qubit* nonnull inttoptr (i64 2 to %Qubit*))
  tail call void @__quantum__qis__cnot__body(%Qubit* nonnull inttoptr (i64 2 to %Qubit*), %Qubit* null)
  tail call void @__quantum__qis__t__adj(%Qubit* nonnull inttoptr (i64 2 to %Qubit*))
  tail call void @__quantum__qis__t__body(%Qubit* null)
  tail call void @__quantum__qis__cnot__body(%Qubit* nonnull inttoptr (i64 1 to %Qubit*), %Qubit* null)
  tail call void @__quantum__qis__h__body(%Qubit* nonnull inttoptr (i64 2 to %Qubit*))
  tail call void @__quantum__qis__x__body(%Qubit* null)
  tail call void @__quantum__qis__h__body(%Qubit* nonnull inttoptr (i64 2 to %Qubit*))
  tail call void @__quantum__qis__x__body(%Qubit* nonnull inttoptr (i64 2 to %Qubit*))
  tail call void @__quantum__qis__h__body(%Qubit* null)
  tail call void @__quantum__qis__x__body(%Qubit* null)
  tail call void @__quantum__qis__z__body(%Qubit* nonnull inttoptr (i64 1 to %Qubit*))
  tail call void @__quantum__qis__cnot__body(%Qubit* null, %Qubit* nonnull inttoptr (i64 1 to %Qubit*))
  tail call void @__quantum__qis__z__body(%Qubit* nonnull inttoptr (i64 1 to %Qubit*))
  tail call void @__quantum__qis__x__body(%Qubit* null)
  tail call void @__quantum__qis__h__body(%Qubit* null)
  tail call void @__quantum__qis__mz__body(%Qubit* null, %Result* null)
  tail call void @__quantum__qis__mz__body(%Qubit* nonnull inttoptr (i64 1 to %Qubit*), %Result* nonnull inttoptr (i64 1 to %Result*))
  ret void
}

declare void @__quantum__qis__cnot__body(%Qubit*, %Qubit*) local_unnamed_addr

declare void @__quantum__qis__z__body(%Qubit*) local_unnamed_addr

declare void @__quantum__qis__t__body(%Qubit*) local_unnamed_addr

declare void @__quantum__qis__t__adj(%Qubit*) local_unnamed_addr

declare void @__quantum__qis__x__body(%Qubit*) local_unnamed_addr

declare void @__quantum__qis__h__body(%Qubit*) local_unnamed_addr

declare void @__quantum__qis__mz__body(%Qubit*, %Result*) local_unnamed_addr

attributes #0 = { "InteropFriendly" "requiredQubits"="3" "requiredResults"="2" }

