```flow
st=>start: N-body snap
op1=>operation: Read snap and separate each components (Host, satellite, disk)
op2=>operation: Compute satellite unbound particles
op3=>operation: Combine host and satellite unbound particles
op4=>operation: Compute BFE coefficients for each component
op5=>operation: Descart noisy coefficients
e=>end: write coefficients

st->op1->op2->op3->op3->op4->op5->e

```

```flow
st=>start: Start
op=>operation: Your Operation
cond=>condition: Yes or No?
e=>end

st->op->cond
cond(yes)->e
cond(no)->op
```
