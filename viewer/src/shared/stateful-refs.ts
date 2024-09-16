// This function exposes a react ref-like interface. However, your components will re-render if this ref is updated

import React from "react";

// https://medium.com/the-non-traditional-developer/creating-a-stateful-ref-object-in-react-fcd56d9dea58
export const useStatefulRef = <Type>(initialVal: Type) => {
    let [currentVal, setCurrentVal] = React.useState<Type>(initialVal);

    const [statefulRef] = React.useState({
        // this getter and setter are just dummy ones
        get current() {
            return initialVal;
        },
        set current(_value: Type) {},
    });
    Object.defineProperty(statefulRef, "current", {
        get: (): Type => currentVal,
        set: (newValue: Type) => {
            if (!Object.is(currentVal, newValue)) {
                currentVal = newValue;
                setCurrentVal(newValue);
            }
        },
    });
    return statefulRef;
};
